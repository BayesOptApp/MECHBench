import os
import shutil
import subprocess
from pathlib import Path
from typing import Union

from src.sob.physical_models.utils.run_openradioss import RunOpenRadioss
from src.sob.physical_models.utils.run_openradioss_docker import (
    RunOpenCourantDocker,
)


def run_OpenRadioss(
    input_file_path: Union[str, Path],
    batch_file_path: Union[str, Path],
    write_vtk: bool = True,
    write_csv: bool = True,
    runStarter: bool = False,
    d3plot: bool = False,
    nt_int: int = 1,
    np_int: int = 1,
):
    """
    Run OpenRadioss locally or OpenCourant through Docker/Podman.

    Parameters
    ----------
    input_file_path:
        Path to the starter input deck, normally <jobname>_0000.rad.

    batch_file_path:
        Use "docker" or "podman" for the container backend.
        Otherwise, provide the local OpenRadioss installation path.

    write_vtk:
        Convert animation files to VTK after the engine runs.

    write_csv:
        Convert time-history files to CSV after the engine runs.

    runStarter:
        If True, run only the starter.

    d3plot:
        Retained for compatibility with the original interface.

    nt_int:
        Number of OpenMP threads.

    np_int:
        Number of MPI processes.

    Notes
    -----
    Existing results associated with the job are deleted before the run,
    as in the original implementation.
    """
    input_file = Path(input_file_path).expanduser().resolve()

    backend = (
        os.fspath(batch_file_path)
        if batch_file_path is not None
        else None
    )

    container_backend = backend in ("docker", "podman")

    command = [
        str(input_file),
        str(nt_int),
        str(np_int),
        "dp",
        "yes" if write_vtk else "no",
        "yes" if write_csv else "no",
        "yes" if runStarter else "no",
        "yes" if d3plot else "no",
    ]

    if container_backend:
        # Pass the application by name. The third positional parameter
        # of RunOpenCourantDocker is the image name.
        run_env = RunOpenCourantDocker(
            command,
            debug=0,
            container_app=backend,
        )

        if shutil.which(backend) is None:
            raise FileNotFoundError(
                f"Cannot find {backend!r} on the Python process PATH. "
                "Check the installation or restart the IDE after "
                "updating PATH."
            )
    else:
        run_env = RunOpenRadioss(
            command,
            0,
            backend,
        )

    if not input_file.is_file():
        raise FileNotFoundError(
            f"Input file not found: {input_file}"
        )

    if not Path(run_env.running_directory).is_dir():
        raise FileNotFoundError(
            f"Run directory not found: {run_env.running_directory}"
        )

    def execute(args, phase):
        """Run a command and report its output if it fails."""


        #print("Command:", repr(args))
        #print("Working directory:", repr(run_env.running_directory))
        #print("Directory exists:", os.path.isdir(run_env.running_directory))

        #print("Docker in runner environment:", shutil.which(
        #    "docker",
        #    path=run_env.environment().get("PATH"),
        #))
        
        if not isinstance(args, (list, tuple)) or not args:
            raise TypeError(
                f"{phase} command must be a nonempty list of arguments"
            )

        try:
            result = subprocess.run(
                args=args,
                env=run_env.environment(),
                cwd=run_env.running_directory,
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Could not start {phase}. "
                f"Executable: {args[0]!r}; "
                f"working directory: {run_env.running_directory!r}. "
                "Check that both exist and the executable is on PATH."
            ) from exc

        if result.returncode != 0:
            output = result.stdout.decode(
                "utf-8",
                errors="replace",
            )

            display_command = subprocess.list2cmdline(
                list(map(str, args))
            )

            raise RuntimeError(
                f"{phase} failed with exit code {result.returncode}.\n"
                f"Command: {display_command}\n"
                f"Output:\n{output}"
            )

        return result

    # Match the original runner's cleanup behavior.
    run_env.delete_previous_results()

    # Run the starter.
    starter_command = run_env.get_starter_command()
    execute(starter_command, "Starter")

    if runStarter:
        return

    # Run the engine decks in order.
    engine_list = run_env.get_engine_input_file_list(repair=True)

    if not engine_list:
        raise FileNotFoundError(
            f"No engine input decks found for {run_env.jobname!r} "
            f"in {run_env.running_directory}"
        )

    for engine_file in engine_list:
        engine_command = run_env.get_engine_command(
            engine_file,
            full_path=True,
        )
        execute(engine_command, f"Engine ({engine_file})")

    # Convert time-history files when requested.
    if write_csv:
        for th_file in run_env.get_th_list():
            run_env.convert_th_to_csv(th_file)

    # Convert animation files when requested.
    if write_vtk:
        for anim_file in run_env.get_animation_list():
            run_env.convert_anim_to_vtk(anim_file)