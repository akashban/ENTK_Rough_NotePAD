import os
import os.path
import sys

def main():
  
# this code will be executed after untarring jump.tar.gz

# check if the simulation has terminated by simply checking if output file exists    

    if os.path.exists('table_1.gro'): # table_1.gro is name of output file

       os.system ('tar -czvf file.tar.gz table_1.gro grompp_B.mdp')
        # make a tar file with the output of stage 1 and a mdp file which has a higher timestep
       print ('lets be agreesive')

    else: 
       os.system ('tar -czvf file.tar.gz conf.gro grompp_A.mdp')
        # make a tar file with initial input file and a mdp file which has a lower timestep

       print ('lets be conservative')



if __name__ == "__main__":
    main()
