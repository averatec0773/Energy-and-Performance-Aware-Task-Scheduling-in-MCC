import matplotlib.pyplot as plt
import os

def read_tasks_from_file(file_name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, file_name)
    
    tasks_data_1_1 = []
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split(", ")
            task = (int(parts[0]), int(parts[1]), int(parts[2]), parts[3].strip('"'))
            tasks_data_1_1.append(task)
    return tasks_data_1_1

def visualize_task_timeline_with_labels(tasks):

    max_time = max(task[2] for task in tasks)
    fig_width = max(16, max_time / 2)
    fig, ax = plt.subplots(6, 1, figsize=(fig_width, 12), sharex=True, constrained_layout=True)
    
    labels = ["Core 1", "Core 2", "Core 3", "Wireless Sending", "Cloud", "Wireless Receiving"]
    colors = ["#FFA07A", "#20B2AA", "#87CEFA", "#FFD700", "#9370DB", "#FF69B4"]
    
    for i in range(6):
        ax[i].set_yticks([])
        ax[i].set_ylabel(labels[i], rotation=0, labelpad=50, fontsize=12, verticalalignment="center")
        ax[i].set_xlim(0, max_time + 5)
        ax[i].grid(True, axis="x", linestyle="--", alpha=0.5)
        ax[i].set_ylim(0, 1)
    
    for task in tasks:
        task_id, ready_time, finish_time, unit = task
        if unit.startswith("Core"):
            core_idx = int(unit.split()[-1]) - 1
            ax[core_idx].barh(0.5, finish_time - ready_time, left=ready_time, height=0.4, 
                              color=colors[core_idx], edgecolor='black')
            ax[core_idx].text(ready_time + (finish_time - ready_time) / 2, 0.5, 
                              f"TASK ID:{task_id}\n{ready_time}-{finish_time}",
                              ha='center', va='center', fontsize=8, color='black')
        elif unit == "Cloud":
            ax[4].barh(0.5, finish_time - ready_time, left=ready_time, height=0.4, 
                       color=colors[4], edgecolor='black')
            ax[4].text(ready_time + (finish_time - ready_time) / 2, 0.5, 
                       f"TASK ID:{task_id}\n{ready_time}-{finish_time}",
                       ha='center', va='center', fontsize=8, color='black')
        elif unit == "Wireless Sending":
            ax[3].barh(0.5, finish_time - ready_time, left=ready_time, height=0.4, 
                       color=colors[3], edgecolor='black')
            ax[3].text(ready_time + (finish_time - ready_time) / 2, 0.5, 
                       f"TASK ID:{task_id}\n{ready_time}-{finish_time}",
                       ha='center', va='center', fontsize=8, color='black')
        elif unit == "Wireless Receiving":
            ax[5].barh(0.5, finish_time - ready_time, left=ready_time, height=0.4, 
                       color=colors[5], edgecolor='black')
            ax[5].text(ready_time + (finish_time - ready_time) / 2, 0.5, 
                       f"TASK ID:{task_id}\n{ready_time}-{finish_time}",
                       ha='center', va='center', fontsize=8, color='black')
    
    ax[0].set_title("Task Execution Timeline with Labels", fontsize=16)
    plt.xlabel("Time", fontsize=14)
    plt.show()

tasks_data_1_1 = read_tasks_from_file("./5_2.txt")
visualize_task_timeline_with_labels(tasks_data_1_1)
