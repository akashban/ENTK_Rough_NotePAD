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

    # Create a Stage object
    s1 = Stage()

    # Create a Task object which creates a file named 'output.txt' of size 1 MB
    t1 = Task()
    os.system('echo START -----------------------------')

 

   
    t1.executable = '/usr/local/gromacs/bin/./gmx grompp -f /home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/grompp.mdp -c /home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/conf.gro -p /home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/topol.top -o /home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/table.tpr -n /home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/index.ndx'
    t1.arguments = ['']

    # Add the Task to the Stage
    s1.add_tasks(t1)

    # Add Stage to the Pipeline
    p.add_stages(s1)

    # Create another Stage object to hold character count tasks
    s2 = Stage()

    # Create a Task object
    t2 = Task()
    t2.executable = '/usr/local/gromacs/bin/./gmx mdrun -deffnm /home/polo-mastaan/Desktop/ENTK/radical.entk/10_15_2019/table -v'
    t2.arguments = ['']
    # Copy data from the task in the first stage to the current task's location ( no point of this line , for now)
    t2.copy_input_data = ['$Pipline_%s_Stage_%s_Task_%s/mdout.mdp' % (p.uid, s1.uid, t1.uid)]

    # Add the Task to the Stage
    s2.add_tasks(t2)

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
