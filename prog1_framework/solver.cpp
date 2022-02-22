#include "solver.h"
using namespace std;

/*************************** DECLARE YOUR HELPER FUNCTIONS HERE ************************/





/*************************** solver.h functions ************************/


void seq_solver(unsigned int n, unsigned int exit_on_first, std::vector<std::vector<unsigned int> >& solns) {

    // Initialize a vector of vector or 2D vector of size 1Xn with value 0
    int rowsSolution = 1;
    std::vector < std::vector<unsigned int> > solutions;

    // TODO: Implement this function
    for (int p1 = 1; p1 <= n; ++p1){

        std::vector<unsigned int> position(n, p1);
        int l = 1;
        int solutionFound = 0;

        // i_th column search from 2_th to n_th columns
        for (int i = 2; i <= n; ++i){   

            // initialize the starting point for position search
            // re-search i_th column position from: position[i-1] + 1
            int startPoint = 1;
            if ((l > n) || (solutionFound == 1)){
                startPoint = position[i-1] + 1;
                solutionFound = 0;
            }

            l = startPoint;

            // position search for i_th column
            for (int k = 0; k < i-1; ++k){
                int j = position[k];

                // threaten strategy
                // if (i == k){
                //     cout << "j="<<j << ", k=" << k << ", i="<<i-1<<", l="<<l <<"\n";
                // }
                if ((j == l) || (std::abs(i - 1 - k) == std::abs(j - l))){
                    k = -1;
                    l += 1;
                    if (l > n){
                        break;
                    }
                }
                // cout << "l=" << l << "  ";
            }

            // for (auto element : position) {
            //         cout << "\n" <<element << ",  ";
            //     } 

            if (l > n){
                i -= 2;
                // cout << "i=" << i+1 << ", ";
            }else{
                position[i-1] = l;
            }
            
            
            if (i == n){
                // for (int m=0; m<n; ++m){
                //     cout<< "solns"<<"\n";
                //     solns[rowsSolution-1][m] = position[m];
                //     cout<< "solns"<<"\n";
                //     cout<< solns[rowsSolution-1][m]<<"\n";
                // }

                cout<< "solns"<<"\n";
                // cout<< position<<"\n";
                // for (int m=0; m<position.size(); ++m){
                //     cout<< position[m]<<"\n";
                //     // position[m] += -1;
                //     int a = position[m];
                //     position[m] = a-1;
                // }

                solns.push_back(position);
                rowsSolution += 1;

                if (exit_on_first == 1){
                    goto end;
                }

                i -= 2;
                solutionFound = 1;
            }
            
            if (i <= 0){
                break;
            }

        }
    }

    

    end:{
        // std::for_each(solns.begin(), solns.end(), [](const auto & elem){elem-1});
        //Print solutions:
        cout << "Number of Solutions=" << rowsSolution-1 << "\n";  

        //solns = &solutions;

        for_each(solns.begin(), solns.end(),
            [](const auto & row ) { 
                for_each(row.begin(), row.end(), 
                        [](const auto & elem){
                            cout << elem << ", ";
                        });
                cout << endl;
            });
        cout << endl;

    }

}


void nqueen_master( unsigned int n,
                    unsigned int k,
                    unsigned int exit_on_first,
                    std::vector<std::vector<unsigned int> >& solns) {
                        
    // TODO: Implement this function

    /* Following is a general high level layout that you can follow
     (you are not obligated to design your solution in this manner.
      This is provided just for your ease). */

    std::vector<unsigned int> final_solution;
    int num_procs;

    // total number of procs in the world
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int rowsSolution = 0;

    int worker = 1;
    int send_times = 0;
    int received_times = 0;
    int tag_id = 1;
    std::vector<unsigned int> recv_final_solution;
    int s = 0;

    for (int p1 = 1; p1 <= n; ++p1){

        vector<int> partial_solution(k, p1);
        int l = 1;
        int solutionFound = 0;

        if (k <= 1){
            partial_solution[0] = p1;
            if (worker <= num_procs-1){
                // cout <<"send_1:" <<"\n"; for_each(partial_solution.begin(), partial_solution.end(),[](const auto & elem){cout <<elem << ", ";}); cout<<" \n" ;
                MPI_Send(&partial_solution[0], k, MPI_INT, worker, tag_id, MPI_COMM_WORLD);
                send_times += 1;
                tag_id += 1;
                solutionFound = 1;
                
            }
            
            worker += 1;

            if (worker > num_procs){
                MPI_Status status;
                int recv_size;
                MPI_Recv(&recv_size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // ----May need waiting time to let worker complete---------
                //MPI_Wait(&req, &status);

                // Specify the rank of responded worker and number of tag:
                int worker_id = status.MPI_SOURCE;
                tag_id = status.MPI_TAG;
                recv_final_solution.resize(recv_size);
                
                MPI_Recv(&recv_final_solution[0], recv_size, MPI_INT, worker_id, tag_id, MPI_COMM_WORLD, &status);
                received_times += 1;

                // cout<<"recv_final_solution"<<"\n"; for_each(recv_final_solution.begin(), recv_final_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n" ;

                // s = recv_final_solution.size()/n;
                // if ((exit_on_first != 1) || ((rowsSolution==1) && (exit_on_first == 1))){
                //     for (int row_i=0; row_i< s; ++row_i){
                //         for (int co_i = 0; co_i < n; ++co_i){
                //             solns[(row_i+rowsSolution)][co_i] = (unsigned int) recv_final_solution[row_i*n + co_i];
                //         }
                //     }

                    
                // }
                // rowsSolution += s;

                s = recv_final_solution.size()/n;
                solns.push_back(recv_final_solution);   
                rowsSolution += s;

                

                // for (int ii = 0; ii < recv_final_solution.size(); ++ii){
                //     solns.push_back(recv_final_solution[ii]);
                // }

                if ((recv_final_solution.size() >= n) && (exit_on_first == 1)){
                    goto exit;
                }

                MPI_Send(&partial_solution[0], k, MPI_INT, worker_id, tag_id, MPI_COMM_WORLD);

                send_times += 1;
                solutionFound = 1;
            }

        }

        if (k != 1){
            // i_th column search from 2_th to k_th columns
            for (int i = 2; i <= k; ++i){   

                // initialize the starting point for position search
                // re-search i_th column postion from: position[i-1] + 1
                int startPoint = 1;
                if ((l > n) || (solutionFound == 1)){
                    startPoint = partial_solution[i-1] + 1;
                    solutionFound = 0;
                }

                l = startPoint;

                // position search for i_th column
                for (int kk = 0; kk < i-1; ++kk){
                    int j = partial_solution[kk];

                    // threaten strategy
                    if ((j == l) || (std::abs(i - 1 - kk) == std::abs(j - l))){
                        kk = -1;
                        l += 1;
                        if (l > n){
                            break;
                        }
                    }
                }

                if (l > n){
                    i -= 2;
                }else{
                    partial_solution[i-1] = l;
                }
                
                // cout<< i <<","<< k<<"," <<worker<<","<<num_procs-1<<"\n";

                if ((i == k) && (worker <= num_procs-1)){
                    /******************* STEP 1: Send one partial solution to each worker ********************/
                    /*
                    * for (all workers) {
                    *      - create a partial solution.
                    *      - send that partial solution to a worker
                    * }
                    */

                    // send that partial solution to a worker

                    // cout <<"send_1:" <<"\n"; for_each(partial_solution.begin(), partial_solution.end(),[](const auto & elem){cout <<elem << ", ";}); cout<<" \n" ;

                    MPI_Send(&partial_solution[0], k, MPI_INT, worker, tag_id, MPI_COMM_WORLD);
                        
                    send_times += 1;
                    tag_id += 1;
                    i -= 2;
                    solutionFound = 1;
                    worker += 1;
                }

                

                if ((i == k) && (worker > num_procs-1)){
                    /******************* STEP 2: Send partial solutions to workers as they respond ********************/
                    /*
                    * while() {
                    *      - receive completed work from a worker processor.
                    *      - create a partial solution
                    *      - send that partial solution to the worker that responded
                    *      - Break when no more partial solutions exist and all workers have responded with jobs handed to them, or if exiting on first solution
                    * }
                    */
                    // 1. receive completed work from any worker processor.
                    // Receives any one of the messages. ''''No order guarante needed''''.
                    // Using Non-blocking receive:
                    MPI_Status status;
                    //MPI_Request req;
                    
                    int recv_size;
                    MPI_Recv(&recv_size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // ----May need waiting time to let worker complete---------
                    //MPI_Wait(&req, &status);

                    // Specify the rank of responded worker and number of tag:
                    int worker_id = status.MPI_SOURCE;
                    tag_id = status.MPI_TAG;
                    recv_final_solution.resize(recv_size);
                    
                    MPI_Recv(&recv_final_solution[0], recv_size, MPI_INT, worker_id, tag_id, MPI_COMM_WORLD, &status);
                    received_times += 1;

                    // cout<<"recv_final_solution"<<"\n"; for_each(recv_final_solution.begin(), recv_final_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n" ;

                    // if (tag_id == 2){
                    //     cout <<"recv_final_solution:" <<"\n"; for_each(recv_final_solution.begin(), recv_final_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n" ;
                    // }

                    // for_each(recv_final_solution.begin(), recv_final_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n" ;

                    // cout << recv_final_solution.size()<<"\n";
                    //store the final solutions

                    // if ((exit_on_first != 1) || ((rowsSolution==1) && (exit_on_first == 1))){
                    //     for (int row_i=0; row_i< s; ++row_i){
                    //         for (int co_i = 0; co_i < n; ++co_i){
                    //             cout << recv_final_solution[row_i*n + co_i]<<"\n";
                    //             cout << row_i+rowsSolution<<"\n";
                    //             cout << co_i<<"\n";
                    //             solns[row_i+rowsSolution][co_i] = 1;
                    //             cout << "solns"<<solns[row_i+rowsSolution][co_i]<<"\n";
                    //             solns[row_i+rowsSolution][co_i] = recv_final_solution[row_i*n + co_i];
                    //             cout << "solns"<<solns[row_i+rowsSolution][co_i]<<"\n";
                    //         }
                    //     }

                        
                    // }
                    

                    // cout << "rowsSolution= "<<rowsSolution<<"\n";
                    // for (int ii = 0; ii < recv_final_solution.size(); ++ii){
                    //     solns.push_back(recv_final_solution[ii]);   
                    // }
                    //cout << "p1= "<<p1<<"\n";
                    s = recv_final_solution.size()/n;
                    solns.push_back(recv_final_solution);   
                    rowsSolution += s;
                    // cout << "rowsSolution= "<<rowsSolution<<"\n";

                    if ((recv_final_solution.size() >= n) && (exit_on_first == 1)){
                        goto exit;
                    }

                    // cout <<"final_solution_update:" <<""; for_each(final_solution.begin(), final_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n" ;

                    // 2 & 3. create a partial solution & send that partial solution to the worker that responded
                    MPI_Send(&partial_solution[0], k, MPI_INT, worker_id, tag_id, MPI_COMM_WORLD);

                    send_times += 1;
                    i -= 2;
                    solutionFound = 1;

                    // cout <<"partial_solution:" <<"\n"; for_each(partial_solution.begin(), partial_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";


                }

                if (i <= 0){
                    break;
                }

            }


        }


    }

    exit:{
        cout<<"send_times="<< send_times<<"\n";
        cout<<"recive_times="<< received_times<<"\n";
            
        }


    while (send_times != received_times){
        // 4. Break when no more partial solutions exist and all workers have responded with jobs handed to them, or if exiting on first solution
        //todo if exiting on first solution---------------
        //当所无任务分配时，查看worker的status， 等待接收最后一波solution
        MPI_Status status;
        //MPI_Request req;
        int recv_size;
        MPI_Recv(&recv_size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // ----May need waiting time to let worker complete---------
        //MPI_Wait(&req, &status);

        // Specify the rank of responded worker and number of tag:
        int worker_id = status.MPI_SOURCE;
        tag_id = status.MPI_TAG;
        recv_final_solution.resize(recv_size);
        MPI_Recv(&recv_final_solution[0], recv_size, MPI_INT, worker_id, tag_id, MPI_COMM_WORLD, &status);
        received_times += 1;
        
        // cout<<"recv_final_solution"<<"\n"; for_each(recv_final_solution.begin(), recv_final_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n" ;

        //store the final solutions
        // s = recv_final_solution.size()/n;
        // if ((exit_on_first != 1) || ((rowsSolution==1) && (exit_on_first == 1))){
        //     for (int row_i=0; row_i< s; ++row_i){
        //         for (int co_i = 0; co_i < n; ++co_i){
        //             solns[row_i+rowsSolution][co_i] = recv_final_solution[row_i*n + co_i];
        //         }
        //     }

            
        // }
        // rowsSolution += s;
        if ((exit_on_first != 1) || ((rowsSolution==1) && (exit_on_first == 1))){
            s = recv_final_solution.size()/n;
            solns.push_back(recv_final_solution);   
            rowsSolution += s;
        }

        // for (int ii = 0; ii < recv_final_solution.size(); ++ii){
        //     final_solution.push_back(recv_final_solution[ii]);
        // }

        // cout <<"fianl solution:" <<"\n"; for_each(final_solution.begin(), final_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";
    }

    

    /********************** STEP 3: Terminate **************************
    *
    * Send a termination/kill signal to all workers.
    *
    */
    std::vector<int> kill(k, -1);
    for (int kill_signal = 1; kill_signal <= num_procs-1; ++kill_signal){
        MPI_Send(&kill[0], k, MPI_INT, kill_signal, kill_signal, MPI_COMM_WORLD);
    }
    goto finalPrint;


    finalPrint: {
        //rowsSolution = final_solution.size() / n;
        cout << "fianl solution #"<< rowsSolution <<"\n";

        for_each(solns.begin(), solns.end(),
            [](const auto & row ) { 
                for_each(row.begin(), row.end(), 
                        [](const auto & elem){
                            cout << elem << ", ";
                        });
                cout << endl;
            });
        cout << endl;


        //std::vector<int> position;
        // for (int row=0; row < rowsSolution; ++row){
        //     cout<<"row="<<row<<"\n";
        //     for (int i = 0; i < n; ++i){
        //         position.push_back(final_solution[n*row+i]);
        //         cout <<"fianl solution:" <<"\n"; for_each(position.begin(), position.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";
        //     }
            // solns[row] = position;
            
            //solns.push_back(position)
            // for (int ii = 0; ii < n; ++ii){
            //     solns.at(row).push_back(position[ii]);
            //     //solns.at(row).at(ii) = position[ii]; 
            //     cout << solns.at(row).at(ii) << ", ";
            // }

        //     if ((row == 0) && (exit_on_first == 1)){
        //         break;
        //     }
        // }

        // cout <<"fianl solution:" <<"\n"; for_each(solns.begin(), solns.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";
        //cout << "send times #"<< send_times<<"\n";
        // cout << "fianl solution #"<< solns.size()/n<<"\n";
        // cout << "fianl solution #"<< final_solution.size()/n<<"\n";
        //cout << "fianl solution 0 = #"<< final_solution<<"\n";
        //cout << "fianl solution 0 = #"<< solns<<"\n";
        //Print solutions:


        // Print final soultion-------------------------------
        // cout <<"fianl solution:"; 
        // for (int i = 0; i < final_solution.size(); ++i){
        //     if (remainder((i), n)==0){
        //         cout << "\n";
        //     }
        //     cout << final_solution[i] << ", ";
        // }

        // solns = final_solution;

        // cout <<"fianl solution:"; 
        // for (int i = 0; i < solns.size(); ++i){
        //     if (remainder((i), n)==0){
        //         cout << "\n";
        //     }
        //     cout << solns[i] << ", ";
        // }

    }
}

void nqueen_worker( unsigned int n,
                    unsigned int k,
                    unsigned int exit_on_first) {
    // TODO: Implement this function

    // Following is a general high level layout that you can follow (you are not obligated to design your solution in this manner. This is provided just for your ease).
    
    int tag_id;
    int master = 0;
    std::vector<int> partial_solution;

    // initialize solution vector
    std::vector<int> solutions;
    std::vector<int> complete_solution(n);
    /*******************************************************************
     *
     * while() {
     *      
     *      wait for a message from master
     *
     *      if (message is a partial job) {
     *              - finish the partial solution
     *              - send all found solutions to master
     *      }
     *
     *      if (message is a kill signal) {
     *
     *              quit
     *
     *      }
     *  }
     */

    while(true) {

        // my proc id
        int proc_id;
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        // receive a message from master
        MPI_Status status;
        partial_solution.resize(k); 
        MPI_Recv(&partial_solution[0], k, MPI_INT, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // cout <<"recive_1:" <<"\n"; for_each(partial_solution.begin(), partial_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n" ;

        // cout << "partial_solution=" << partial_solution << "\n";
        tag_id = status.MPI_TAG;

        if (partial_solution[0] == -1){
            break;
        }

        
        // if (message is a partial job)
        if (partial_solution[0] != -1){

            //- finish the partial solution
            // initialize
            int l = 1;
            int solutionFound = 0;
            
            for (int ii = 0; ii <= k-1; ++ii){
                complete_solution[ii] = partial_solution[ii];
            }

            // i_th column search from 5_th to n_th columns
            for (int i = k+1; i <= n; ++i){   

                // initialize the starting point for position search
                // re-search i_th column postion from: position[i-1] + 1
                int startPoint = 1;
                if ((l > n) || (solutionFound == 1)){
                    startPoint = complete_solution[i-1] + 1;
                    solutionFound = 0;
                }

                // position search for i_th column
                l = startPoint;
                for (int kk = 0; kk < i-1; ++kk){
                    int j = complete_solution[kk];

                    // threaten strategy
                    if ((j == l) || (std::abs(i - 1 - kk) == std::abs(j - l))){
                        kk = -1;
                        l += 1;
                        if (l > n){
                            break;
                        }
                    }
                }

                if (l > n){
                    i -= 2;
                }else{
                    complete_solution[i-1] = l;
                }
                
                // find one solution, then store in 'solutions' vector, continum search:
                if (i == n){

                    // cout <<"complete solutions:" <<"\n"; for_each(complete_solution.begin(), complete_solution.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";

                    for (int ii = 0; ii < n; ++ii){
                        solutions.push_back(complete_solution[ii]);

                        // cout <<"complete solutions:" <<"\n"; for_each(solutions.begin(), solutions.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";
                    }

                    solutionFound = 1;
                    i -= 2;   
                }
                
                // find all solutions & send the 'solutions' vertor to 'master'
                // tags changes to proc_id, because master need to finish first step first, and than receive completed work from a worker processor
                // if (tag_id == 1){
                //     cout <<"i=" <<i<<"\n";
                // }
                
                if (k == 1){
                    if (i < k){
                        int send_size = solutions.size();

                        // if (tag_id == 1){
                        //     cout << "solution_size="<< solutions.size()<<"\n";
                        // }

                        // if (tag_id == 1){
                        // cout <<"send complete solutions:" <<"\n"; for_each(solutions.begin(), solutions.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";}

                        // cout <<"send complete solutions:" <<"\n"; for_each(solutions.begin(), solutions.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";
                        
                        MPI_Send(&send_size, 1, MPI_INT, master, tag_id, MPI_COMM_WORLD);
                        MPI_Send(&solutions[0], send_size, MPI_INT, master, tag_id, MPI_COMM_WORLD);
                        solutions.clear();
                        break;
                    }
                }

                if (i < k-1){
                    int send_size = solutions.size();

                    // if (tag_id == 1){
                    //     cout << "solution_size="<< solutions.size()<<"\n";
                    // }

                    // if (tag_id == 1){
                    // cout <<"send complete solutions:" <<"\n"; for_each(solutions.begin(), solutions.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";}

                    // cout <<"send complete solutions:" <<"\n"; for_each(solutions.begin(), solutions.end(),[](const auto & elem){cout << elem << ", ";}); cout<<" \n";
                    

                    MPI_Send(&send_size, 1, MPI_INT, master, tag_id, MPI_COMM_WORLD);
                    MPI_Send(&solutions[0], send_size, MPI_INT, master, tag_id, MPI_COMM_WORLD);
                    solutions.clear();
                    break;
                }

            }
        }

    }
}



/*************************** DEFINE YOUR HELPER FUNCTIONS HERE ************************/


