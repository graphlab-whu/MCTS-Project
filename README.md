# MCTS-Project

## ABSTRACT

The recently proposed historical $k$-core query introduces a new paradigm of structure analysis for temporal graphs.
However, the query processing based on the existing PHC-index, which preserves the distinct ''core time'' of each vertex, needs to traverse all vertices for each query, even though the results usually contain only a small subset of vertices.
Inspired by the traditional $k$-shell that ensures the optimal $k$-core query processing, we propose a novel concept called ''core time shell'', which reveals the hierarchical structure of vertices with respect to their core time.
Based on the core time shell, we design a time-space balanced Merged Core Time Shell index (MCTS-index).
It is theoretically guaranteed that, the MCTS-index provides the approximately optimal query performance, and has the approximately same space complexity as the PHC-index.
Moreover, we leverage the MCTS-index to efficiently address the brand-new ''when'' historical $k$-core queries orthogonal to the current ''what'' historical $k$-core queries.
Our experimental results on ten real-world temporal graphs demonstrate both the superior efficiency of processing ''what'' queries and the effectiveness of processing versatile ''when'' queries for the MCTS-index.

## File Description
1. phc.cpp: Source code of PHC-index
2. mcts(id).cpp: Source code of the (id) version, which ranks the vertices in each core time shell by their ids
3. mcts(lp).cpp: Source code of the (lp) version, which heuristically moves a vertex to the last position of its new core time shell
4. mcts.cpp: Constructs the MCTS-index and generates test data to conduct experiments and compare query efficiency with PHC-index and TCD-Single.
## Documentation

Experiments are conducted on a Linux machine with a 2.2 GHz CPU and 120 GB of RAM. At least 32GB memory is recommended. All algorithms are implemented using the C++ Standard Template Library and compiled with the g++ compiler at the -O3 optimization level. The confirmed steps of compiling the source code are provided below:

1. Configure the g++ compilation environment in Linux.
2. Add the source code file and compile it to generate the executable file.
3. Add the graph file to the directory of executable file.
4. Run the executable file with command.

The command to compile the source file and generate an executable file is as follows:
The command to construct the PHC-index is:  ```.\[phc exe file] [graph file] [phc file]```, where ```[phc exe file]``` is the executable file generated from phc.cpp, [graph file] is the file storing the graph, and the constructed PHC-index will be stored in ```[phc file]```.

The command to construct the MCTS-index(id) is:  ```.\[mcts(id) exe file] [phc file]```, where ```[mcts(id) exe file]``` is the executable file generated from mcts(id).cpp, and [phc file] is the file storing the PHC-index.

The command to construct the MCTS-index(lp) is:  ```.\[mcts(lp) exe file] [phc file]```, where ```[mcts(lp) exe file]``` is the executable file generated from mcts(lp).cpp, and [phc file] is the file storing the PHC-index.

The command to construct the MCTS-index and conduct the experiment is: ```.\[mcts exe file] [graph file] [phc file]```, where ```[mcts exe file]``` is the executable file generated from mcts.cpp, [graph file] is the file storing the graph, and [phc file] is the file storing the PHC-index.

The dataset used in the experiment can be downloaded from Google Drive:
**link**: https://drive.google.com/drive/folders/1jcD4jIYIprxq-gqL6vuG8zXinb8-i68H

The format of the dataset storing the graph is as follows: each line contains three integers u, v, and t, representing a connection between node u and node v at time t.

Since constructing the PHC-index is time-consuming, we provide a pre-constructed PHC-index for convenience.
The link is shared below:
**link**: https://pan.baidu.com/s/1G0dyH4JMqMiW2WlzCENd2w?pwd=3yci

**password**: 3yci

The storage format of the PHC-index is as follows: the first line contains three numbers vern, kmax, and tmax, representing the number of nodes, the maximum coreness, and the maximum timestamp, respectively.
The second line contains vern integers, representing the IDs of each node in the original graph.
Next, for each of the vern nodes, the first line stores the node's maximum coreness in the graph. For each of the possible k values up to the maximum coreness, the first line stores the number t, which indicates the count of pairs ```[ts, te]``` corresponding to k in the PHC-index. The following t lines each contain two numbers, representing a pair ```[ts, te]```.


## Contact
If you have any questions, contact us by sending an email to clock@whu.edu.cn / zhiwang@whu.edu.cn
