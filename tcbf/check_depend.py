import os
import shutil
import sys
from tcbf.run_command import run_command

def Check_dependencies():
    print("Checking dependency! ")
    dependency_program = "mash minimap2 mcl lastdb".split()
    have_command_not_in_path = False
    for command in dependency_program:
        if shutil.which(command) is None:
            if os.path.exists(os.path.join("external", command)):
                os.environ["PATH"] = os.environ["PATH"] + ":" + os.path.join(os.path.abspath("."),"external")
                os.environ["PATH"] = os.environ["PATH"] + ":" + os.path.join(os.path.abspath("."), "external","bin")
                continue
            print(f"{command} is not in PATH !!!")
            have_command_not_in_path = True
            is_downloading = input(f"auto download {command} to ./external?  Y/N")
            if is_downloading.upper() == "Y" or is_downloading.upper() == "YES":
                print(f"安装{command}中")
                download_dependency(command)

            else:
                sys.exit()
    if have_command_not_in_path:
        os.environ["PATH"] = os.environ["PATH"] + ":" + os.path.join(os.path.abspath("."),"external")






def download_minimap2():
    command = "wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 -O minimap2.tar.bz2;" \
          "tar xvjf  minimap2.tar.bz2;" \
          "cp minimap2-2.24_x64-linux/minimap2 external"
    run_command(command)


def download_mcl():

    command = "wget https://raw.githubusercontent.com/davidemms/OrthoFinder/master/scripts_of/bin/mcl;" \
              "chmod +x mcl;" \
              "mv mcl external"
    run_command(command)


def download_mash():
    command = "wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar;" \
              "tar xf mash-Linux64-v2.3.tar;" \
              "cp mash-Linux64-v2.3/mash external"
    run_command(command)

def download_last():
    install_dir = os.path.join(os.path.abspath("."),"external")
    command = "git clone https://gitlab.com/mcfrith/last.git;" \
              "cd last;" \
              "make;" \
              f"make install prefix={install_dir}"

    run_command(command)
    os.environ["PATH"] = os.environ["PATH"] + ":" + os.path.join(os.path.abspath("."), "external", "bin")

def download_gffread():
    install_dir = os.path.join(os.path.abspath("."), "external")
    command = "https://github.com/gpertea/gffread.git" \
              "cd gffread;" \
              "make release;" \
              f"mv gffread {install_dir}"
    run_command(command)
def download_dependency(command):
    if not os.path.exists("external"):
        os.mkdir("external")
    if command == "minimap2":
        download_minimap2()
    elif command == "mcl":
        download_mcl()
    elif command == "mash":
        download_mash()
    elif command == "lastdb":
        download_last()

