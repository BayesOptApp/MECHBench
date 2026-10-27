# Copyright 1986-2025 Altair Engineering Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
# IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# ====================================================================================
# This is a modification performed by Ivan Olarte-Rodriguez 
# 24/04/2025 <-> 14:05
# Since the idea is top operate OpenRadioss from afar, then the OpenRadioss path is 
# altered such that the user should just point to the main path of the OpenRadioss
# environment and run the case afterwards.
# ====================================================================================
"""Run OpenCourant through Docker or Podman using the native runner interface.

Uses the official convenience commands documented at:
https://opencourant.org/downloads/
Starter/engine use convenience names; converters use the explicit Linux x86-64
names supplied by the user. The converter names target an x86-64 image.
Requires Python 3.9+ and Docker/Podman running Linux containers.
"""

import os
import re
import glob
import shlex
import shutil
import subprocess
from pathlib import Path, PurePosixPath
from typing import Optional, Sequence, Union


def get_deck_runid(file):
    match = re.fullmatch(r".+_([0-9]{4})\.rad", Path(file).name, re.IGNORECASE)
    if Path(file).suffix.lower() == ".rad" and match is None:
        raise ValueError("RAD decks must be named <jobname>_<four digits>.rad")
    return int(match.group(1)) if match else 0


class RunOpenCourantDocker:
    """Use command=[deck, nt, np, precision, ..., ..., starter_only].

    container_app selects docker or podman. Optionally set
    container_executable to the absolute path to the chosen executable.
    """

    def __init__(self, command: Sequence, debug: int,
                 image: str = "opencourant/opencourant",
                 mount_point: str = "/work", container_app: str = "docker",
                 container_executable: Optional[Union[str, Path]] = None):
        if len(command) < 7:
            raise ValueError("command must contain at least seven entries")
        source = Path(command[0]).expanduser().resolve()
        self.extension = source.suffix.lower()
        if self.extension == ".rad":
            self.run_id = get_deck_runid(source)
            self.jobname = source.stem[:-5]
        elif self.extension == ".inp":
            self.jobname, self.run_id = source.stem, 0
            # Preserve the native runner's converted-deck convention.
            source = source.with_name(self.jobname + "_0000.rad")
        elif self.extension == ".k":
            self.jobname, self.run_id = source.stem, 0
            # Preserve the native runner's converted-deck convention.
            source = source.with_name(self.jobname + ".k")
        else:
            raise ValueError("Input must be a .rad, .k or .inp deck")
        self._file = str(source.as_posix())
        self.working_dir = source.parent
        self.running_directory = str(self.working_dir.as_posix())
        self.nt = str(self._positive_int(command[1], "nt"))
        self.np = str(self._positive_int(command[2], "np"))
        self.precision = str(command[3]).lower()
        if self.precision not in ("sp", "dp"):
            raise ValueError("precision must be 'sp' or 'dp'")
        flag = command[6]
        if isinstance(flag, str):
            flag = flag.strip().lower()
            if flag not in ("true", "false", "1", "0", "yes", "no", ""):
                raise ValueError("starter_only must be a boolean or boolean string")
            flag = flag in ("true", "1", "yes")
        self.starter_only = bool(flag)
        self.debug = debug
        self.image = image
        mount = PurePosixPath(mount_point)
        if not mount.is_absolute() or ".." in mount.parts:
            raise ValueError("mount_point must be an absolute container path")
        self.mount_point = str(mount)
        self.container_app = container_app
        self.container_executable = container_executable
        self.environment()

    @staticmethod
    def _positive_int(value, name):
        result = int(value)
        if result < 1 or str(value).strip() != str(result):
            raise ValueError(f"{name} must be a positive integer")
        return result

    @property
    def nb_threads(self):
        return int(self.nt)

    @property
    def nb_mpi(self):
        return int(self.np)

    @property
    def file(self):
        return self._file

    @property
    def container_app(self):
        return self._container_app

    @container_app.setter
    def container_app(self, value):
        if value not in ("docker", "podman"):
            raise ValueError("container_app must be 'docker' or 'podman'")
        self._container_app = value

    def environment(self):
        """Return the host launch environment; image owns solver library paths."""
        self.custom_env = os.environ.copy()
        #self.custom_env.update(OMP_NUM_THREADS=self.nt, OMP_STACKSIZE="400m",
        #                       KMP_STACKSIZE="400m")
        return self.custom_env

    def _resolve_container_executable(self):
        """Return an absolute host executable path for Docker or Podman.

        An explicit container_executable overrides PATH discovery. It must be
        an absolute path to a file. Otherwise, the selected application is
        discovered using the current Python process PATH.
        """
        if self.container_executable is not None:
            executable = Path(self.container_executable).expanduser()
            if not executable.is_absolute() or not executable.is_file():
                raise FileNotFoundError(
                    f"Container executable must be an existing absolute file path: {executable}"
                )
            return str(executable.resolve())
        executable = shutil.which(self.container_app)
        if executable is None:
            raise FileNotFoundError(
                f"Cannot find {self.container_app!r} on PATH. Install it, "
                "restart the Python process, or supply container_executable."
            )
        return str(Path(executable).resolve())

    def _check_runner(self):
        try:
            executable = self._resolve_container_executable()
            return subprocess.run([executable, "info"],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                env=self.custom_env).returncode == 0
        except OSError:
            return False

    def pull_image(self):
        executable = self._resolve_container_executable()
        subprocess.run([executable, "pull", self.image], check=True,
                       env=self.custom_env)

    def _container_file(self, file, full_path=False):
        path = Path(file)
        # Relative names are relative to the mounted run directory.
        host = path.resolve() if path.is_absolute() else (self.working_dir / path).resolve()
        try:
            relative = host.relative_to(self.working_dir)
        except ValueError as exc:
            raise ValueError("File must be inside the mounted run directory") from exc
        target = PurePosixPath(*relative.parts)
        return str(PurePosixPath(self.mount_point) / target) if full_path else str(target)

    def _docker_cmd(self, command: Union[str, Sequence[str]], solver_env=True):
        # Keep the official entrypoint and pass argument lists without a shell.
        # Shell snippets are used only for explicit executable availability checks.
        arguments = ["sh", "-c", command] if isinstance(command, str) else list(command)
        cmd = [self._resolve_container_executable(), "run", "--rm", "-v",
               f"{self.working_dir.as_posix()}:{self.mount_point}",
               "-w", self.mount_point]
        if solver_env:
            cmd += ["-e", f"OMP_NUM_THREADS={self.nt}"]
        return cmd + [self.image] + arguments

    def _print_command(self, command):
        if self.debug:
            # Windows quoting is for display only; subprocess receives a list.
            print(subprocess.list2cmdline(command) if os.name == "nt"
                  else shlex.join(command))

    def _run(self, command, log_file: Optional[str] = None):
        cmd = self._docker_cmd(command)
        self._print_command(cmd)
        if log_file:
            with (self.working_dir / log_file).open("w") as stream:
                return subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT,
                                      env=self.custom_env).returncode
        return subprocess.run(cmd, env=self.custom_env).returncode

    def _executable(self, kind, precision=None):
        precision = self.precision if precision is None else precision
        mpi = "_ompi" if kind == "engine" and self.nb_mpi > 1 else ""
        suffix = "_sp" if precision == "sp" else ""
        return f"{kind}{mpi}{suffix}"

    def get_starter_command(self, full_path=False):
        return self._docker_cmd([self._executable("starter"), "-i",
            self._container_file(self.file, full_path), "-np", self.np])

    def get_engine_command(self, engine_input, full_path=False):
        command = [self._executable("engine"), "-i",
                   self._container_file(engine_input, full_path)]
        if self.nb_mpi > 1:
            command = ["mpiexec", "-np", self.np] + command
        return self._docker_cmd(command)

    def run_starter(self, rootname=None, input_suffix="_0000.rad"):
        deck = self.file if rootname is None else rootname + input_suffix
        return self._run([self._executable("starter"), "-i",
            self._container_file(deck), "-np", self.np],
            f"{self.jobname}_starter.log")

    def run_engine(self, rootname=None, run_number="0001"):
        rootname = self.jobname if rootname is None else rootname
        number = f"{int(run_number):04d}"
        command = self.get_engine_command(f"{rootname}_{number}.rad")
        self._print_command(command)
        with (self.working_dir / f"{rootname}_{number}_engine.log").open("w") as stream:
            return subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                  env=self.custom_env).returncode

    def run(self, rootname=None, starter_suffix="_0000.rad", run_number=None):
        """Run the starter and all engine decks; stop on the first failure.

        A nonzero input deck starts directly with that engine deck. An explicit
        run_number selects only that engine. .inp inputs require an already
        converted <jobname>_0000.rad, as in the native runner.
        """
        if not self._check_runner():
            raise RuntimeError("Container application or daemon is unavailable")
        root = self.jobname if rootname is None else rootname
        if self.run_id == 0 or rootname is not None:
            deck = self.file if rootname is None else str(self.working_dir / (root + starter_suffix))
            if not Path(deck).is_file():
                raise FileNotFoundError(deck)
            if self.run_starter(rootname, starter_suffix) != 0:
                return False
        if self.starter_only:
            return True
        if run_number is not None:
            return self.run_engine(root, run_number) == 0
        pattern = re.compile(re.escape(root) + r"_([0-9]{4})\.rad")
        decks = sorted((int(m.group(1)), p) for p in self.working_dir.iterdir()
                       if p.is_file() and (m := pattern.fullmatch(p.name))
                       and int(m.group(1)) >= max(1, self.run_id))
        if not decks:
            raise FileNotFoundError(f"No engine decks found for {root}")
        return all(self.run_engine(root, number) == 0 for number, _ in decks)

    def list_results(self, rootname=None):
        root = self.jobname if rootname is None else rootname
        return sorted(p for p in self.working_dir.iterdir()
                      if p.is_file() and p.name.startswith(root)
                      and p.suffix.lower() != ".rad")

    def get_jobname_runid_rundirectory(self):
        return self.jobname, self.run_id, self.running_directory

    def get_engine_input_file_list(self, repair=False):
        pattern = re.compile(re.escape(self.jobname) + r"_([0-9]{4})\.rad")
        decks = sorted(p for p in self.working_dir.iterdir() if p.is_file()
                       and (m := pattern.fullmatch(p.name)) and int(m.group(1)) > 0)
        if repair:
            for deck in decks:
                text = deck.read_text()
                cards = {line.strip() for line in text.splitlines()}
                missing = [c for c in ("/TH/TITLE", "/ANIM/ELEM/THIC") if c not in cards]
                if missing:
                    # Insert before /END so the solver reads the added cards.
                    lines = text.splitlines(keepends=True)
                    index = next((i for i, line in enumerate(lines)
                                  if line.strip() == "/END"), len(lines))
                    if index and not lines[index - 1].endswith("\n"):
                        lines[index - 1] += "\n"
                    lines[index:index] = [card + "\n" for card in missing]
                    deck.write_text("".join(lines))
        return [p.name for p in decks]

    def get_animation_list(self):
        pattern = re.compile(re.escape(self.jobname) + r"A[0-9]{3,}")
        return sorted(p.name for p in self.working_dir.iterdir()
                      if p.is_file() and pattern.fullmatch(p.name))

    def get_th_list(self):
        pattern = re.compile(re.escape(self.jobname) + r"T[0-9]{2}")
        return sorted(p.name for p in self.working_dir.iterdir()
                      if p.is_file() and pattern.fullmatch(p.name))

    def get_anim_to_vtk_command(self, anim_file, full_path=False):
        """Return the command; redirect its stdout to the desired VTK file."""
        return self._docker_cmd(["anim_to_vtk_linux64_gf",
            self._container_file(anim_file, full_path)], solver_env=False)

    def get_th_to_csv_command(self, th_file, full_path=False):
        return self._docker_cmd(["th_to_csv_linux64_gf",
            self._container_file(th_file, full_path)], solver_env=False)

    def convert_anim_to_vtk(self, anim_file):
        """Capture converter stdout as <anim_file>.vtk in the run folder."""
        relative = self._container_file(anim_file)
        output = self.working_dir / (relative + ".vtk")
        temporary = output.with_name(output.name + ".tmp")
        command = self.get_anim_to_vtk_command(anim_file)
        self._print_command(command)
        try:
            # Binary capture preserves exact output bytes on Windows too.
            with temporary.open("wb") as stream:
                proc = subprocess.run(command, stdout=stream,
                                      stderr=subprocess.PIPE, env=self.custom_env)
            if proc.returncode:
                raise RuntimeError("Animation conversion failed: " +
                                   (proc.stderr or b"" ).decode(errors="replace"))
            temporary.replace(output)
        finally:
            temporary.unlink(missing_ok=True)
        return output

    def convert_th_to_csv(self, th_file):
        relative = self._container_file(th_file)
        command = self.get_th_to_csv_command(th_file)
        self._print_command(command)
        proc = subprocess.run(command,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=self.custom_env)
        if proc.returncode:
            raise RuntimeError("Time-history conversion failed: " + (proc.stdout or ""))
        output = self.working_dir / (relative + ".csv")
        if not output.is_file():
            raise RuntimeError(f"Output file not created: {output}")
        return output

    def _check_precision(self, precision):
        names = [self._executable("starter", precision), self._executable("engine", precision)]
        script = " && ".join("command -v " + shlex.quote(name) + " >/dev/null" for name in names)
        return self._run(script) == 0

    def check_sp_exes(self):
        return self._check_precision("sp")

    def check_dp_exes(self):
        return self._check_precision("dp")

    def delete_previous_results(self):
       
        extensions_to_delete = [".h3d", "_[0-9][0-9][0-9][0-9].out", "_[0-9][0-9][0-9][0-9].ctl", "_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].rst",
                                "T[0-9][0-9]", "T[0-9][0-9].csv", "A[0-9][0-9][0-9]", "A[0-9][0-9][0-9][0-9]",
                                "A[0-9][0-9][0-9].vtk", "A[0-9][0-9][0-9][0-9].vtk", ".d3plot", ".d3plot[0-9][0-9]",
                                ".d3plot[0-9][0-9][0-9]", ".d3plot[0-9][0-9][0-9][0-9]",'.tmp']

        # Delete files with specified extensions
        for ext_pattern in extensions_to_delete:
            files_to_delete = os.path.join(self.running_directory, self.jobname + f"{ext_pattern}")
            matching_files = [file for file in glob.glob(files_to_delete) if re.fullmatch(re.escape(self.jobname) + ext_pattern.replace('.', r'\.'), os.path.basename(file))]
            #print(matching_files)
            for file in matching_files:
                os.remove(file)

