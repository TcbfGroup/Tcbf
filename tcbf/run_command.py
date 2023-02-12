import subprocess
from multiprocessing import Pool
def run_command(command):
    """
    此段代码借鉴了  109-130行  https://github.com/davidemms/OrthoFinder/blob/master/scripts_of/__main__.py
    command : 需要执行的命令
    """

    capture = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    stdout, stderr = capture.communicate()
    try:
        stdout = stdout.decode()
        stderr = stderr.decode()
    except(UnicodeError, AttributeError):
        stdout = stdout.encode()
        stderr = stderr.encode()
    n_stderr_lines = stderr.count('\n')
    if capture.returncode != 0:
        print(f"Returned error :{capture.returncode}")
        print(f"Command: {command}")
        print("stdout:\n-------")
        print(stdout)
        if n_stderr_lines > 0:
            print("stderr:\n-------")
            print(stderr)
        raise Exception("error")
    return stdout


def parall_run(commands,process_number):
    if process_number <= 0:
        from multiprocessing import cpu_count
        process_number = cpu_count()
    with Pool(process_number)as p:
       p.map(run_command,commands)
