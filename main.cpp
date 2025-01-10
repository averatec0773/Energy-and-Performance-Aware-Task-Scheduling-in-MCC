// PROJECT 2
// YIYANG HUANG

#include <iostream>
#include <cfloat>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <stack>
#include <fstream>
#include <string>

using namespace std;

// Global Value
const int Core_num = 3;
const double T_max_factor = 1.5;
const double T_s = 3;
const double T_c = 1;
const double T_r = 1;
const double T_remote = T_s + T_c + T_r;

// Task
struct Task {
    double id;                      // Task ID
    double localTime[Core_num];     // Execution time on the 3 local cores
    double priority = -1;           // Priority
    double remote = -1;             // -1: Unassigned, 0: Local, 1: Uploaded to the cloud
    double assignedUnit = -1;       // -1: Unassigned, 0~2: Core 1~3, Core_num: Uploaded to the cloud
    double global_finish_time = -1; // Global finish time
    double readyTimeLocal = -1;     // Ready time for local task
    double finishTimeLocal = -1;    // Finish time for local task
    double readyTimeRemote[3] = {-1, -1, -1};  // Ready times for cloud tasks (upload, execute, receive)
    double finishTimeRemote[3] = {-1, -1, -1}; // Finish times for cloud tasks (upload, execute, receive)
};

// Unit
struct Unit {
    double id = -1;                // 0-Core_num-1: Local cores, Core_num: Cloud core
    double e = -1;                 // Energy consumption
    double taskid = -1;            // Currently executing task ID
    double last_finishTime = 0;    // Last finish time
    double last_s = 0;             // Last start time
    double last_c = 0;             // Last computation time
    double last_r = 0;             // Last receive time
    vector<int> arrange_taskid;    // List of arranged task IDs
};

// Print Task Info
void print_task(const vector<Task>& tasks, double& T_total, double& E_total) {
    cout << "T_total: " << T_total << endl;
    cout << "E_total: " << E_total << endl;
    cout << "Task Information:\n";
    cout << "Task ID | Priority | Assign State | Assigned Unit | Ready Time | Finish Time\n";
    cout << "------------------------------------------------------------\n";

    for (const auto& task : tasks) {
        cout << "Task " << task.id << " | "
             << task.priority << " | ";

        if (task.remote == 0){
            cout << "Local" << " | ";
        } else {
            cout << "Cloud" << " | ";
        }
        if (task.assignedUnit == Core_num) {
            cout << "Cloud | "
                 << "R: " << task.readyTimeRemote[0] << "/" << task.readyTimeRemote[1] << "/" << task.readyTimeRemote[2] << " | "
                 << "F: " << task.finishTimeRemote[0] << "/" << task.finishTimeRemote[1] << "/" << task.finishTimeRemote[2];
        } else {
            cout << "Core " << (task.assignedUnit + 1) << " | "
                 << "R: " << task.readyTimeLocal << " | "
                 << "F: " << task.finishTimeLocal;
        }
        cout << "\n";
    }
    cout << "------------------------------------------------------------\n";
}

// Save Task Info into Local File
void saveTasksToFile(const vector<Task>& tasks, const string &filename) {
    ofstream outFile(filename, ios::out | ios::trunc);
    if (!outFile) {
        cerr << "Error: Unable to open or create file!" << endl;
        return;
    }

    for (int i = 0; i < tasks.size(); ++i) {
        const Task &task = tasks[i];
        if (task.assignedUnit >= 0 && task.assignedUnit < Core_num) {
            outFile << task.id << ", " << task.readyTimeLocal << ", "
                    << task.finishTimeLocal << ", \"Core " << task.assignedUnit + 1 << "\"\n";
        } else if (task.assignedUnit == Core_num) {
            const char *labels[] = {"Wireless Sending", "Cloud", "Wireless Receiving"};
            for (int j = 0; j < 3; ++j) {
                outFile << task.id << ", " << task.readyTimeRemote[j] << ", "
                        << task.finishTimeRemote[j] << ", \"" << labels[j] << "\"\n";
            }
        }
    }

    outFile.close();
    cout << "Tasks saved to " << filename << " successfully!" << endl;
}

// Step One: Initial Scheduling
// Phase 1：Primary Assignment
void primary_assignment(vector<Task>& tasks) {
    for (auto& task : tasks) {
        double local_time = min({task.localTime[0], task.localTime[1], task.localTime[2]});
        double cloud_time = T_remote;

        if (cloud_time < local_time) {
            task.remote = 1;
        } else {
            task.remote = 0;
        }
    }
}

// Phase 2：Task Prioritizing
double calculate_priority(int taskId, vector<Task>& tasks, const vector<vector<int>>& adjMatrix, vector<bool>& visited) {
    if (visited[taskId]) {
        return tasks[taskId].priority;
    }

    double w_i = 0;
    if (tasks[taskId].assignedUnit == Core_num) {
        w_i = T_remote;
    } else {
        w_i = (tasks[taskId].localTime[0] + tasks[taskId].localTime[1] + tasks[taskId].localTime[2]) / Core_num;
    }

    bool isExitTask = true;
    double maxChildPriority = 0;
    for (int i = 0; i < tasks.size(); ++i) {
        if (adjMatrix[taskId][i] == 1) {
            isExitTask = false;
            maxChildPriority = max(maxChildPriority, calculate_priority(i, tasks, adjMatrix, visited));
        }
    }

    if (isExitTask) {
        tasks[taskId].priority = w_i;
    } else {
        tasks[taskId].priority = w_i + maxChildPriority;
    }

    visited[taskId] = true;
    return tasks[taskId].priority;
}

vector<int> task_prioritizing(vector<Task>& tasks, const vector<vector<int>>& adjMatrix) {
    vector<bool> visited(tasks.size(), false);
    for (int i = 0; i < tasks.size(); ++i) {
        calculate_priority(i, tasks, adjMatrix, visited);
    }

    vector<pair<double, int>> task_priority_pairs;
    for (const auto& task : tasks) {
        task_priority_pairs.push_back({task.priority, task.id});
    }

    sort(task_priority_pairs.begin(), task_priority_pairs.end(),
         [](const pair<double, int>& a, const pair<double, int>& b) {
             return a.first > b.first;
         });

    vector<int> priority_queue;
    for (const auto& pair : task_priority_pairs) {
        priority_queue.push_back(pair.second);
    }

    // for (int i = 0; i < priority_queue.size(); ++i)
    //     cout << priority_queue[i] << " ";
    // cout << endl;

    return priority_queue;
}

// Phase 3：Execution Unit Selection
void execution_unit_selection(const vector<int>& priority_queue, vector<Task>& tasks, vector<Unit>& units, vector<vector<int>>& adjMatrix, double& T_min, double& E_total){
    T_min = 0;

    int task_index = priority_queue[0] - 1;
    tasks[task_index].readyTimeLocal = 0;
    tasks[task_index].readyTimeRemote[0] = 0;
    
    // Remote
    if (tasks[task_index].remote == 1) {
        tasks[task_index].finishTimeRemote[0] = tasks[task_index].readyTimeRemote[0] + T_s;
        tasks[task_index].readyTimeRemote[1] = tasks[task_index].finishTimeRemote[0];
        tasks[task_index].finishTimeRemote[1] = tasks[task_index].readyTimeRemote[1] + T_c;
        tasks[task_index].readyTimeRemote[2] = tasks[task_index].finishTimeRemote[1];
        tasks[task_index].finishTimeRemote[2] = tasks[task_index].readyTimeRemote[2] + T_r;

        tasks[task_index].finishTimeLocal = 0;
        tasks[task_index].global_finish_time = tasks[task_index].finishTimeRemote[2];
        units[Core_num].last_s = tasks[task_index].finishTimeRemote[0];
        units[Core_num].last_c = tasks[task_index].finishTimeRemote[1];
        units[Core_num].last_r = tasks[task_index].finishTimeRemote[2];
        units[Core_num].arrange_taskid.push_back(task_index);
        tasks[task_index].assignedUnit = Core_num;

        E_total += T_s * units[Core_num].e;
    } 
    // Local
    else {
        double local_time_min = DBL_MAX;
        int core_index;
        for (int i = 0; i < Core_num; i++) {
            if (tasks[task_index].localTime[i] < local_time_min) {
                local_time_min = tasks[task_index].localTime[i];
                core_index = i;
            }
        }

        double Local_estimate = tasks[task_index].readyTimeLocal + local_time_min;
        double Remote_estimate = tasks[task_index].readyTimeRemote[0] + T_s + T_c + T_r;

        // Keep Local
        if (Local_estimate <= Remote_estimate) {
            tasks[task_index].finishTimeLocal = Local_estimate;
            tasks[task_index].global_finish_time = tasks[task_index].finishTimeLocal;
            units[core_index].last_finishTime = tasks[task_index].global_finish_time;
            units[core_index].arrange_taskid.push_back(task_index);
            tasks[task_index].assignedUnit = core_index;

            E_total += tasks[task_index].localTime[core_index] * units[core_index].e;
        } 
        // Change Remote
        else {
            tasks[task_index].finishTimeRemote[0] = tasks[task_index].readyTimeRemote[0] + T_s;
            tasks[task_index].readyTimeRemote[1] = tasks[task_index].finishTimeRemote[0];
            tasks[task_index].finishTimeRemote[1] = tasks[task_index].readyTimeRemote[1] + T_c;
            tasks[task_index].readyTimeRemote[2] = tasks[task_index].finishTimeRemote[1];
            tasks[task_index].finishTimeRemote[2] = Remote_estimate;

            tasks[task_index].finishTimeLocal = 0;
            tasks[task_index].global_finish_time = tasks[task_index].finishTimeRemote[2];
            units[Core_num].last_s = tasks[task_index].finishTimeRemote[0];
            units[Core_num].last_c = tasks[task_index].finishTimeRemote[1];
            units[Core_num].last_r = tasks[task_index].finishTimeRemote[2];
            units[Core_num].arrange_taskid.push_back(task_index);
            tasks[task_index].assignedUnit = Core_num;
            tasks[task_index].remote = 1;

            E_total += T_s * units[Core_num].e;
        }
    }

    for(int p = 1; p < tasks.size(); ++p){
        int task_index = priority_queue[p] - 1;

        // Max Local Ready
        double pre_local_time = 0;
        for(int j=0; j<tasks.size(); j++){
            if(adjMatrix[j][task_index] == 1 && pre_local_time < max(tasks[j].finishTimeLocal,tasks[j].finishTimeRemote[2])){
                pre_local_time = max(tasks[j].finishTimeLocal,tasks[j].finishTimeRemote[2]);
            }
        }
        tasks[task_index].readyTimeLocal = pre_local_time;

        // Max Sending Ready
        double pre_sending_time = 0;
        for(int j=0; j<tasks.size(); j++){
            if(adjMatrix[j][task_index] == 1 && pre_sending_time < max(tasks[j].finishTimeLocal,tasks[j].finishTimeRemote[0])){
                pre_sending_time = max(tasks[j].finishTimeLocal,tasks[j].finishTimeRemote[2]);
            }
        }
        tasks[task_index].readyTimeRemote[0] = max(pre_sending_time,units[Core_num].last_s);
        
        // Max Computing Ready
        double pre_computing_time = 0;
        for(int j=0; j<tasks.size(); j++){
            if(adjMatrix[j][task_index] == 1 && pre_computing_time < tasks[j].finishTimeRemote[1]){
                pre_computing_time = tasks[j].finishTimeRemote[2] + T_s;
            }
        }
        tasks[task_index].readyTimeRemote[1] = max(pre_computing_time,units[Core_num].last_c);

        // Remote
        if(tasks[task_index].remote == 1){
            tasks[task_index].finishTimeRemote[0] = tasks[task_index].readyTimeRemote[0] + T_s;
            tasks[task_index].readyTimeRemote[1] = tasks[task_index].finishTimeRemote[0];
            tasks[task_index].finishTimeRemote[1] = tasks[task_index].readyTimeRemote[1] + T_c;
            tasks[task_index].readyTimeRemote[2] = tasks[task_index].finishTimeRemote[1];
            tasks[task_index].finishTimeRemote[2] = tasks[task_index].readyTimeRemote[2] + T_r;

            tasks[task_index].finishTimeLocal = 0;
            tasks[task_index].global_finish_time = tasks[task_index].finishTimeRemote[2];
            units[Core_num].last_s = tasks[task_index].finishTimeRemote[0];
            units[Core_num].last_c = tasks[task_index].finishTimeRemote[1];
            units[Core_num].last_r = tasks[task_index].finishTimeRemote[2];
            units[Core_num].arrange_taskid.push_back(task_index);
            tasks[task_index].assignedUnit = Core_num;

            E_total += T_s * units[Core_num].e;
        }
        // Local
        else{
            int ready_time, core_index;
            double new_local_time = DBL_MAX;
            for(int j=0; j < Core_num; j++){
                ready_time = max(tasks[task_index].readyTimeLocal,units[j].last_finishTime);
                if(new_local_time > ready_time + tasks[task_index].localTime[j]){
                    new_local_time = ready_time + tasks[task_index].localTime[j];
                    core_index = j;
                }
            }
            tasks[task_index].readyTimeLocal = new_local_time - tasks[task_index].localTime[core_index];

            double Local_estimate = new_local_time;
            double Remote_estimate = tasks[task_index].readyTimeRemote[0] + T_s + T_c +  T_r;

            // Keep Local
            if(Local_estimate <= Remote_estimate){
                tasks[task_index].finishTimeLocal = Local_estimate;
                tasks[task_index].global_finish_time = tasks[task_index].finishTimeLocal;
                units[core_index].last_finishTime = tasks[task_index].global_finish_time;
                units[core_index].arrange_taskid.push_back(task_index);
                tasks[task_index].assignedUnit = core_index;

                E_total += tasks[task_index].localTime[core_index] * units[core_index].e;
            }
            // Change Remote
            else{
                // cout << "Task " << task_index + 1 << "Change Remote" << endl;
                tasks[task_index].finishTimeRemote[0] = tasks[task_index].readyTimeRemote[0] + T_s;
                tasks[task_index].readyTimeRemote[1] = tasks[task_index].finishTimeRemote[0];
                tasks[task_index].finishTimeRemote[1] = tasks[task_index].readyTimeRemote[1] + T_c;
                tasks[task_index].readyTimeRemote[2] = tasks[task_index].finishTimeRemote[1];
                tasks[task_index].finishTimeRemote[2] = Remote_estimate;

                tasks[task_index].finishTimeLocal = 0;
                tasks[task_index].global_finish_time = tasks[task_index].finishTimeRemote[2];
                units[Core_num].last_s = tasks[task_index].finishTimeRemote[0];
                units[Core_num].last_c = tasks[task_index].finishTimeRemote[1];
                units[Core_num].last_r = tasks[task_index].finishTimeRemote[2];
                units[Core_num].arrange_taskid.push_back(task_index);
                tasks[task_index].assignedUnit = Core_num;
                tasks[task_index].remote = 1;

                E_total += T_s * units[Core_num].e;
            }
        }

        T_min = max(T_min, tasks[task_index].global_finish_time);

    }
}

vector<Task> initial_scheduling(vector<Task>& tasks, vector<Unit>& units, vector<vector<int>>& adjMatrix, vector<vector<int>>& Seq, double& T_min, double& E_total) {
    primary_assignment(tasks);
    vector<int> priority_queue = task_prioritizing(tasks, adjMatrix);
    execution_unit_selection(priority_queue, tasks, units, adjMatrix, T_min, E_total);

    for (size_t core_index = 0; core_index < units.size(); ++core_index) {
        Seq[core_index] = units[core_index].arrange_taskid;
    }

    return tasks;
}

// Step Two: Task Migration
vector<vector<int>> create_seq(vector<Task>& tasks, vector<vector<int>>& adjMatrix, vector<vector<int>>& Seq, int old_core, int new_core, int& task_id){

    vector<vector<int>> tmp_Seq = Seq;
    auto &oldQueue = tmp_Seq[old_core];

    if (find(oldQueue.begin(), oldQueue.end(), task_id) == oldQueue.end()) {
        
        cerr << "Error: Task " << task_id << " not found in Core " << old_core << " queue for removal." << endl;
        cerr << "Current queue: ";
        for (const auto& task : oldQueue) {
            cerr << task << " ";
        }
        cerr << endl;
        throw runtime_error("Task not found for removal in create_seq.");
    }

    oldQueue.erase(remove(oldQueue.begin(), oldQueue.end(), task_id), oldQueue.end());

    double task_start_time = tasks[task_id].readyTimeLocal;
    auto &newQueue = tmp_Seq[new_core];
    bool inserted = false;

    for (auto it = newQueue.begin(); it != newQueue.end(); ++it) {
        int current_task_id = *it;

        double current_task_readytime = max(tasks[current_task_id].readyTimeRemote[0],tasks[current_task_id].readyTimeLocal);

        if (current_task_readytime > task_start_time) {
            newQueue.insert(it, task_id);
            inserted = true;
            break;
        }
        else if (current_task_readytime == task_start_time && adjMatrix[current_task_id][task_id] != 1){
            newQueue.insert(it+1, task_id);
            inserted = true;
            break;
        }
        else if (current_task_readytime == task_start_time && adjMatrix[current_task_id][task_id] == 1){
            newQueue.insert(it, task_id);
            inserted = true;
            break;
        }
    }

    if (!inserted) {
        newQueue.push_back(task_id);
    }

    int original_tasks = 0;
    for (const auto& queue : Seq) {
        original_tasks += queue.size();
    }
    int total_tasks = 0;
    for (const auto& queue : tmp_Seq) {
        total_tasks += queue.size();
    }

    return tmp_Seq;

}

vector<Task> kernel(vector<Task>& tasks, vector<Unit>& units, vector<vector<int>>& adjMatrix, vector<vector<int>>& Seq, double& current_T, double& current_E){
    vector<int> ready1(tasks.size(), 0);
    vector<int> ready2(tasks.size(), 0);
    vector<bool> scheduled(tasks.size(), false);

    // Initialization ready1
    for (int i=0;i<tasks.size();++i){
        for (int j=0;j<tasks.size();++j){
            if (adjMatrix[i][j]==1){
                ready1[j] += 1;
            }
        }
    }

    // Initialization ready2
    for (int i=0;i<Seq.size();++i){
        for (int j=0;j<Seq[i].size();++j){
            ready2[Seq[i][j]] += j;
        }
    }

    stack<int> taskStack;
    for (int i = 0; i < tasks.size(); ++i) {
        if (ready1[i] == 0 && ready2[i] == 0) {
            tasks[i].readyTimeLocal = 0;
            tasks[i].readyTimeRemote[0] = 0;
            taskStack.push(i);
            scheduled[i] = true;
            // cout << i << " ";
        }
    }

    while (!taskStack.empty()) {

        int task_index = taskStack.top();
        taskStack.pop();
        int assign_core = tasks[task_index].assignedUnit;

        if (tasks[task_index].assignedUnit != Core_num) {

            // Max Local Ready
            double pre_local_time = 0;
            for(int j=0; j<tasks.size(); j++){
                if(adjMatrix[j][task_index] == 1 && pre_local_time < max(tasks[j].finishTimeLocal,tasks[j].finishTimeRemote[2])){
                    pre_local_time = max(tasks[j].finishTimeLocal,tasks[j].finishTimeRemote[2]);
                }
            }
            tasks[task_index].readyTimeLocal = max(pre_local_time,units[assign_core].last_finishTime);

            double local_time = tasks[task_index].localTime[assign_core];

            tasks[task_index].finishTimeLocal = tasks[task_index].readyTimeLocal + local_time;
            tasks[task_index].global_finish_time = tasks[task_index].finishTimeLocal;
            units[assign_core].last_finishTime = tasks[task_index].global_finish_time;

            current_E += tasks[task_index].localTime[assign_core] * units[assign_core].e;
            // cout<<"TASK: " << task_index+1 << "CURRENT E: " << current_E << endl;

        } else {
            // Max Sending Ready
            double pre_sending_time = 0;
            for(int j=0; j<tasks.size(); j++){
                if(adjMatrix[j][task_index] == 1 && pre_sending_time < max(tasks[j].finishTimeLocal,tasks[j].finishTimeRemote[0])){
                    pre_sending_time = max(tasks[j].finishTimeLocal,tasks[j].finishTimeRemote[2]);
                }
            }
            tasks[task_index].readyTimeRemote[0] = max(pre_sending_time,units[Core_num].last_s);

            // Max Computing Ready
            double pre_computing_time = 0;
            for(int j=0; j<tasks.size(); j++){
                if(adjMatrix[j][task_index] == 1 && pre_computing_time < tasks[j].finishTimeRemote[1]){
                    pre_computing_time = tasks[j].finishTimeRemote[2] + T_s;
                }
            }
            tasks[task_index].readyTimeRemote[1] = max(pre_computing_time,units[Core_num].last_c);

            tasks[task_index].finishTimeLocal = 0;
            tasks[task_index].finishTimeRemote[0] = tasks[task_index].readyTimeRemote[0] + T_s;
            tasks[task_index].readyTimeRemote[1] = tasks[task_index].finishTimeRemote[0];
            tasks[task_index].finishTimeRemote[1] = tasks[task_index].readyTimeRemote[1] + T_c;
            tasks[task_index].readyTimeRemote[2] = tasks[task_index].finishTimeRemote[1];
            tasks[task_index].finishTimeRemote[2] = tasks[task_index].readyTimeRemote[2] + T_r;

            tasks[task_index].global_finish_time = tasks[task_index].finishTimeRemote[2];
            units[Core_num].last_s = tasks[task_index].finishTimeRemote[0];
            units[Core_num].last_c = tasks[task_index].finishTimeRemote[1];
            units[Core_num].last_r = tasks[task_index].finishTimeRemote[2];
            tasks[task_index].remote = 1;

            current_E += T_s * units[Core_num].e;
            // cout << "TASK: " << task_index+1 << "CURRENT E: " << current_E << endl;

        }

        current_T = max(current_T, tasks[task_index].global_finish_time);

        for (int i = 0; i < tasks.size(); ++i) {
            if (adjMatrix[task_index][i] == 1) {
                ready1[i] -= 1;
            }
        }

        int index_task_in_seq = -1;
        
        for (int j = 0; j < (int)Seq[assign_core].size(); ++j) {
            if (Seq[assign_core][j] == task_index) {
                index_task_in_seq = j;
                break;
            }
        }

        if (index_task_in_seq != -1) {
            for (int j = index_task_in_seq + 1; j < (int)Seq[assign_core].size(); ++j) {
                int next_task_id = Seq[assign_core][j];
                ready2[next_task_id] -= 1;
            }
        }

        for (int t = 0; t < tasks.size(); t++) {
            if (!scheduled[t] && ready1[t] == 0 && ready2[t] == 0) {
                taskStack.push(t);
                scheduled[t] = true;
                // cout << t << " ";
            }
        }
    }

    return tasks;

}

vector<Task> outer_loop(vector<Task>& tasks, vector<Unit>& units, vector<vector<int>>& adjMatrix, vector<vector<int>>& Seq, double& T_total, double& E_total, double T_max_factor = 1.5) 
{
    vector<Task> best_current_tasks = tasks;
    vector<Unit> best_current_units = units;
    vector<vector<int>> best_current_Seq = Seq;
    double best_current_T = T_total;
    double best_current_E = E_total;

    bool improved = true;
    int w_iter = 0;

    while (improved) {
        if(w_iter == 100){break;}
        improved = false;
        w_iter += 1;

        vector<Task> best_no_increase_tasks = best_current_tasks;
        vector<Unit> best_no_increase_units = best_current_units;
        vector<vector<int>> best_no_increase_Seq = best_current_Seq;
        double best_no_increase_T = best_current_T;
        double best_no_increase_E = best_current_E;
        bool found_no_increase = false;

        vector<Task> best_ratio_tasks = best_current_tasks;
        vector<Unit> best_ratio_units = best_current_units;
        vector<vector<int>> best_ratio_Seq = best_current_Seq;
        double best_ratio_T = best_current_T;
        double best_ratio_E = best_current_E;
        bool found_increase = false;
        
        for (int i = 0; i < (int)best_current_tasks.size(); ++i) {
            int current_unit = best_current_tasks[i].assignedUnit;
            if (best_current_tasks[i].remote == 0 && current_unit < Core_num) {
                // Change Core
                for (int new_core = 0; new_core <= Core_num; ++new_core) {
                    if (new_core == current_unit) continue;

                    // cout << "Try Task " << i+1 << " From Core " << current_unit+1 << " to Core " << new_core+1 << endl; 

                    // Deep Copy
                    vector<Task> temp_tasks = best_current_tasks;
                    vector<Unit> temp_units = best_current_units;
                    vector<vector<int>> temp_Seq = best_current_Seq;
                    double current_T = 0, current_E = 0;

                    for (int core_index = 0; core_index < (int)units.size(); ++core_index) {
                        temp_units[core_index].last_finishTime = 0;
                        temp_units[core_index].last_s = 0;
                        temp_units[core_index].last_c = 0;
                        temp_units[core_index].last_r = 0;
                        temp_units[core_index].arrange_taskid.clear();
                    }

                    int old_core = current_unit;
                    temp_tasks[i].assignedUnit = new_core;

                    if(new_core<Core_num){temp_tasks[i].remote = 0;}
                    else{temp_tasks[i].remote = 0;}

                    vector<vector<int>> current_Seq = create_seq(temp_tasks, adjMatrix, temp_Seq, old_core, new_core, i);

                    // cout<<endl<<"Iter: "<<w_iter<<" Task "<<i+1<<" to core "<<new_core+1<<endl;
                    // print_task(temp_tasks,current_T,current_E);
                    // cout<<"current_Seq: "<<endl;
                    // for(int i=0;i<current_Seq.size();++i){
                    //     cout<<"S"<<"["<<i+1<<"]: ";
                    //     for (int j=0;j<(int)current_Seq[i].size();++j){
                    //         cout << current_Seq[i][j]+1 << " ";
                    //     }
                    //     cout<<endl;
                    // }
                    // cout<<endl;

                    vector<Task> new_solution = kernel(temp_tasks, temp_units, adjMatrix, current_Seq, current_T, current_E);

                    if (current_T <= best_current_T) {
                        if (current_E < best_no_increase_E) {
                            best_no_increase_tasks = new_solution;
                            best_no_increase_units = temp_units;
                            best_no_increase_Seq = current_Seq;
                            best_no_increase_T = current_T;
                            best_no_increase_E = current_E;
                            found_no_increase = true;
                        }
                    } else if (current_T <= T_total * T_max_factor) {
                        if (current_E < best_ratio_E) {
                            best_ratio_tasks = new_solution;
                            best_ratio_units = temp_units;
                            best_ratio_Seq = current_Seq;
                            best_ratio_T = current_T;
                            best_ratio_E = current_E;
                            found_increase = true;
                        }
                    }
                }
            }
        }

        if (found_no_increase) {
            best_current_tasks = best_no_increase_tasks;
            best_current_units = best_no_increase_units;
            best_current_Seq = best_no_increase_Seq;
            best_current_T = best_no_increase_T;
            best_current_E = best_no_increase_E;
            improved = true;
        } else if (found_increase) {
            best_current_tasks = best_ratio_tasks;
            best_current_units = best_ratio_units;
            best_current_Seq = best_ratio_Seq;
            best_current_T = best_ratio_T;
            best_current_E = best_ratio_E;
            improved = true;
        } 

        cout<<endl<<"Best Solution at "<<w_iter<<endl;
        print_task(best_current_tasks,best_current_T,best_current_E);

        // cout << "Current Best Seq state:\n";
        // for (size_t core_index = 0; core_index < best_current_Seq.size(); ++core_index) {
        //     cout << "Core " << core_index + 1 << ": ";
        //     for (const auto& task : best_current_Seq[core_index]) {
        //         cout << "Task " << task + 1 << " ";
        //     }
        //     cout << endl;
        // }

    }

    tasks = best_current_tasks;
    units = best_current_units;
    T_total = best_current_T;
    E_total = best_current_E;
    return best_current_tasks;
}

int main() {

    vector<Task> tasks_1 = {
        {1,  {9, 7, 5}},
        {2,  {8, 6, 5}},
        {3,  {6, 5, 4}},
        {4,  {7, 5, 3}},
        {5,  {5, 4, 2}},
        {6,  {7, 6, 4}},
        {7,  {8, 5, 3}},
        {8,  {6, 4, 2}},
        {9,  {5, 3, 2}},
        {10, {7, 4, 2}}
    };

    vector<Task> tasks_2 = {
        {1,  {9, 7, 5}},
        {2,  {8, 6, 5}},
        {3,  {6, 5, 4}},
        {4,  {7, 5, 3}},
        {5,  {5, 4, 2}},
        {6,  {7, 6, 4}},
        {7,  {8, 5, 3}},
        {8,  {6, 4, 2}},
        {9,  {5, 3, 2}},
        {10, {7, 4, 2}},
        {11, {9, 6, 4}},
        {12, {8, 7, 5}},
        {13, {7, 5, 3}},
        {14, {6, 4, 2}},
        {15, {8, 6, 4}},
        {16, {7, 6, 5}},
        {17, {9, 8, 6}},
        {18, {6, 5, 3}},
        {19, {5, 3, 1}},
        {20, {7, 5, 4}}
    };

    vector<Unit> units = {
        {1, 1},
        {2, 2},
        {3, 4},
        {4, 0.5}
    };

    vector<vector<int>> adjMatrix_10_1, adjMatrix_10_2, adjMatrix_20_1, adjMatrix_20_2, adjMatrix_20_3;

    adjMatrix_10_1 = vector<vector<int>>(tasks_1.size(), vector<int>(tasks_1.size(), 0));
    //                   1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    adjMatrix_10_1[0] = {0, 1, 1, 1, 1, 1, 0, 0, 0, 0}; 
    adjMatrix_10_1[1] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 0}; 
    adjMatrix_10_1[2] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}; 
    adjMatrix_10_1[3] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 0}; 
    adjMatrix_10_1[4] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}; 
    adjMatrix_10_1[5] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}; 
    adjMatrix_10_1[6] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; 
    adjMatrix_10_1[7] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; 
    adjMatrix_10_1[8] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; 
    adjMatrix_10_1[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
    //                   1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    adjMatrix_10_2 = vector<vector<int>>(tasks_1.size(), vector<int>(tasks_1.size(), 0));
    //                   1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    adjMatrix_10_2[0] = {0, 1, 1, 1, 1, 1, 0, 0, 0, 0}; 
    adjMatrix_10_2[1] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}; 
    adjMatrix_10_2[2] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}; 
    adjMatrix_10_2[3] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 0}; 
    adjMatrix_10_2[4] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}; 
    adjMatrix_10_2[5] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}; 
    adjMatrix_10_2[6] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; 
    adjMatrix_10_2[7] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; 
    adjMatrix_10_2[8] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; 
    adjMatrix_10_1[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 

    adjMatrix_20_1 = vector<vector<int>>(tasks_2.size(), vector<int>(tasks_2.size(), 0));
    //                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19,20
    adjMatrix_20_1[0]  = {0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[1]  = {0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[2]  = {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[3]  = {0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[4]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[5]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[6]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[7]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[8]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_1[9]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
    adjMatrix_20_1[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    adjMatrix_20_1[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0};
    adjMatrix_20_1[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    adjMatrix_20_1[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0};
    adjMatrix_20_1[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_1[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_1[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_1[17] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_1[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_1[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    adjMatrix_20_2 = vector<vector<int>>(tasks_2.size(), vector<int>(tasks_2.size(), 0));
    //                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19,20
    adjMatrix_20_2[0]  = {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[1]  = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[2]  = {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[3]  = {0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[4]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[5]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[6]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[7]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[8]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_2[9]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
    adjMatrix_20_2[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    adjMatrix_20_2[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0};
    adjMatrix_20_2[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    adjMatrix_20_2[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0};
    adjMatrix_20_2[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_2[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_2[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_2[17] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_2[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    adjMatrix_20_2[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    adjMatrix_20_3 = vector<vector<int>>(tasks_2.size(), vector<int>(tasks_2.size(), 0));
    //                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19,20
    adjMatrix_20_3[0]  = {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[1]  = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[2]  = {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[3]  = {0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[4]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[5]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[6]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[7]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[8]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[9]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};
    adjMatrix_20_3[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
    adjMatrix_20_3[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0};
    adjMatrix_20_3[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
    adjMatrix_20_3[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    adjMatrix_20_3[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
    adjMatrix_20_3[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
    adjMatrix_20_3[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
    adjMatrix_20_3[17] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    adjMatrix_20_3[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    double T_total = -1;
    double E_total = 0;

    vector<Task> tasks = tasks_2;
    vector<vector<int>> adjMatrix = adjMatrix_20_3;

    vector<vector<int>> Seq(units.size());

    // Step One
    vector<Task> tasks_step1 = initial_scheduling(tasks, units, adjMatrix, Seq, T_total, E_total);
    double T_step1 = T_total;
    double E_step1 = E_total;

    cout << "Step One Completed.\n";
    print_task(tasks_step1, T_step1, E_step1);
    saveTasksToFile(tasks_step1,"5_1.txt");

    // Step Two
    vector<Task> tasks_step2 = outer_loop(tasks, units, adjMatrix, Seq, T_total, E_total, T_max_factor);
    double T_step2 = T_total;
    double E_step2 = E_total;

    cout << "Step Two Completed.\n";
    print_task(tasks_step2, T_step2, E_step2);
    saveTasksToFile(tasks_step2,"5_2.txt");

    return 0;
}
