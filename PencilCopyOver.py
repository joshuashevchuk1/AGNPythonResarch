import os
import sys

def run():
    root = os.getcwd() # root dir is fucked up due to a space
    dir_run_list = next(os.walk('.'))[1]
    maxOrbits = 500
    ivar = 0
    while (ivar <= maxOrbits):
        moveImages(dir_run_list,ivar,root) # copy over images into PencilAnalysis
        ivar = ivar + 50
    os.chdir(root)

def moveImages(dir_run_list,ivar,root):
     for i in range(len(dir_run_list)):
         os.chdir(dir_run_list[i])
         name = os.path.split(os.getcwd())[1]
         fullName = str(name) + "-ivar-" + str(ivar) + ".png"
         os.chdir(root)
         os.system('cp -rf ' + dir_run_list[i] + "/" + str(fullName) + " ./" + "Pencil_Analysis")
         os.chdir(root)

run()  
