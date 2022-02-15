#include "solver.h"
using namespace std;

/*************************** DECLARE YOUR HELPER FUNCTIONS HERE ************************/





/*************************** solver.h functions ************************/


void seq_solver(unsigned int n, unsigned int exit_on_first, std::vector<std::vector<unsigned int> >& solns) {

    // Initialize a vector of vector or 2D vector of size 1Xn with value 0
    int columnsSolution = n;
    int rowsSolution = 1;
    std::vector < vector<int> > solutions;

    // TODO: Implement this function
    for (int p1 = 8; p1 <= 8; ++p1){

        std::vector<int> position(n, p1);
        int l = 1;
        int solutionFound = 0;

        // i_th column search from 2_th to 8_th columns
        for (int i = 2; i <= n; ++i){   

            // initialize the starting point for position search
            // re-search i_th column postion from: position[i-1] + 1
            int startPoint = 1;
            if ((l > 8) ||(solutionFound = 1)){
                startPoint = position[i-1] + 1;
            }

            l = startPoint;

            // position search for i_th column
            for (int k = 0; k < i-1; ++k){
                int j = position[k];

                // threaten strategy
                // if (i == 4){
                //     cout << "j="<<j << ", k=" << k << ", i="<<i-1<<", l="<<l <<"\n";
                // }
                if ((j == l) || (std::abs(i - 1 - k) == std::abs(j - l))){
                    k = -1;
                    l += 1;
                    if (l > 8){
                        break;
                    }
                }
                // cout << "l=" << l << "  ";
            }

            // for (auto element : position) {
            //         cout << "\n" <<element << ",  ";
            //     } 

            if (l > 8){
                i -= 2;
                cout << "i=" << i+1 << ", ";
            }else{
                position[i-1] = l;
            }
            
            // Print solutions:
            if (i == n){
                solutions.push_back(position);
                rowsSolution += 1;
                cout << "rowsSolution=" << rowsSolution << ", ";

                i -= 2;
                solutionFound = 1;

                for_each(solutions.begin(), solutions.end(),
                    [](const auto & row ) { 
                        for_each(row.begin(), row.end(), 
                                [](const auto & elem){
                                    cout << elem << ", ";
                                });
                        cout << endl;
                    });
                cout << endl; 
            }

            if (i <= 0){
                break;
            }
                
        }
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


    /******************* STEP 1: Send one partial solution to each worker ********************/
    /*
     * for (all workers) {
     *      - create a partial solution.
     *      - send that partial solution to a worker
     * }
     */


    /******************* STEP 2: Send partial solutions to workers as they respond ********************/
    /*
     * while() {
     *      - receive completed work from a worker processor.
     *      - create a partial solution
     *      - send that partial solution to the worker that responded
     *      - Break when no more partial solutions exist and all workers have responded with jobs handed to them, or if exiting on first solution
     * }
     */

    /********************** STEP 3: Terminate **************************
     *
     * Send a termination/kill signal to all workers.
     *
     */





}

void nqueen_worker( unsigned int n,
                    unsigned int k,
                    unsigned int exit_on_first) {



    // TODO: Implement this function

    // Following is a general high level layout that you can follow (you are not obligated to design your solution in this manner. This is provided just for your ease).

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


}



/*************************** DEFINE YOUR HELPER FUNCTIONS HERE ************************/







