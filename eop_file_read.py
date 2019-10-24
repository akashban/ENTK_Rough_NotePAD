from radical.entk import Pipeline, Stage, Task, AppManager
import os

# ------------------------------------------------------------------------------
# Set default verbosity
if os.environ.get('RADICAL_ENTK_VERBOSE') == None:
    os.environ['RADICAL_ENTK_REPORT'] = 'True'


# Description of how the RabbitMQ process is accessible
# No need to change/set any variables if you installed RabbitMQ has a system
# process. If you are running RabbitMQ under a docker container or another
# VM, set "RMQ_HOSTNAME" and "RMQ_PORT" in the session where you are running
# this script.
hostname = os.environ.get('RMQ_HOSTNAME', 'localhost')
port = os.environ.get('RMQ_PORT', 5672)


def generate_pipeline():

    # Create a Pipeline object
    p = Pipeline()
    p.name ='p1'

    # Create a Stage object
    s1 = Stage()
    s1.name ='s1'

    # Create a Task1 in stage1
    t1 = Task()
    t1.name ='t1'
    os.system('echo START -----------------------------')
    
    t1.upload_input_data = ['/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/grompp.mdp',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/conf.gro',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/topol.top',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/index.ndx', 
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/table_A_H.xvg',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/table_H_H.xvg',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/table.xvg' ]


    t1.executable = '/usr/local/gromacs/bin/./gmx'

    t1.pre_exec =['/usr/local/gromacs/bin/./gmx' + ' grompp -f grompp.mdp -c conf.gro -p topol.top -o table.tpr -n index.ndx']

   
    t1.arguments = ['mdrun',
                    '-deffnm', 
                    'table',
                    '-v']
    t1.post_exec =[' mv ' + ' table.gro ' + ' table_1.gro ']

    # Add the Task1 to the Stage1
    s1.add_tasks(t1)

    # Create a Task2 in stage1
    t2 = Task()
    t2.name ='t2'

    
    t2.upload_input_data = ['/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/grompp.mdp',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/conf.gro',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/topol.top',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/index.ndx', 
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/table_A_H.xvg',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/table_H_H.xvg',
                            '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/table.xvg' ]


    t2.executable = '/usr/local/gromacs/bin/./gmx'

    t2.pre_exec =['/usr/local/gromacs/bin/./gmx' + ' grompp -f grompp.mdp -c conf.gro -p topol.top -o table.tpr -n index.ndx']

   
    t2.arguments = ['mdrun',
                    '-deffnm', 
                    'table',
                    '-v']
    t2.post_exec =[' mv ' + ' table.gro ' + ' table_1.gro ', ' pwd ']

    # Add the Task2 to the Stage1
    s1.add_tasks(t2)

    # Add Stage to the Pipeline
    p.add_stages(s1)

    # Create another Stage object 
    s2 = Stage()
    s2.name ='s2'

    # Create a Task object in stage2
    t3 = Task()
    t3.name ='t3'
    t3.copy_input_data = ['$Pipline_%s_Stage_%s_Task_%s/grompp.mdp' % (p.name, s1.name, t1.name),
                          '$Pipline_%s_Stage_%s_Task_%s/table_1.gro' % (p.name, s1.name, t1.name) ,
                          '$Pipline_%s_Stage_%s_Task_%s/topol.top' % (p.name, s1.name, t1.name),
                          '$Pipline_%s_Stage_%s_Task_%s/index.ndx' % (p.name, s1.name, t1.name),
                          '$Pipline_%s_Stage_%s_Task_%s/conf.gro' % (p.name, s1.name, t1.name),
                          '$Pipline_%s_Stage_%s_Task_%s/table.xvg' % (p.name, s1.name, t1.name),
                          '$Pipline_%s_Stage_%s_Task_%s/table_A_H.xvg' % (p.name, s1.name, t1.name),
                          '$Pipline_%s_Stage_%s_Task_%s/table_H_H.xvg' % (p.name, s1.name, t1.name)]


    t3.executable = '/usr/local/gromacs/bin/./gmx'

    t3.pre_exec =['/usr/local/gromacs/bin/./gmx' + ' grompp -f grompp.mdp -c table_1.gro -p topol.top -o table.tpr -n index.ndx']

   
    t3.arguments = ['mdrun',
                    '-deffnm', 
                    'table',
                    '-v']
    t3.post_exec =[' mv ' + ' table.gro ' + ' table_2.gro ',
                   ' pwd ']
   

    # Add the Task to the Stage2
    s2.add_tasks(t3)

     # Create  Task2 object in stage2
    t4 = Task()
    t4.name ='t4'
    t4.copy_input_data = ['/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/something.cpp',
                          '/home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/read.py',
                           '$Pipline_%s_Stage_%s_Task_%s/conf.gro' % (p.name, s1.name, t1.name)]


    t4.executable = 'g++ -std=c++11'

  
   
    t4.arguments = ['something.cpp']
    t4.post_exec =['./a.out ', ' python ' + ' read.py ']

 


    
   

    # Add the Task to the Stage2
    s2.add_tasks(t4)




    # Add Stage to the Pipeline
    p.add_stages(s2)


    return p


if __name__ == '__main__':

    pipelines = []

    for cnt in range(1):
        pipelines.append(generate_pipeline())

    # Create Application Manager
    appman = AppManager(hostname=hostname, port=port)

    # Create a dictionary describe four mandatory keys:
    # resource, walltime, and cpus
    # resource is 'local.localhost' to execute locally
    res_dict = {

        'resource': 'local.localhost',
        'walltime': 10,
        'cpus': 1
    }

    # Assign resource request description to the Application Manager
    appman.resource_desc = res_dict

    # Assign the workflow as a set or list of Pipelines to the Application Manager
    # Note: The list order is not guaranteed to be preserved
    appman.workflow = set(pipelines)

    # Run the Application Manager
    appman.run()
