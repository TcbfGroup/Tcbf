import subprocess
import sys
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
        sys.exit()
    return stdout
