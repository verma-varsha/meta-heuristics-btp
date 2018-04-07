#include <bits/stdc++.h>
#include <random>
#include <time.h>
using namespace std;

class Particle{
    public:
        vector<double> pos;
        int makespan;
};

class MachineData{
    public:
        int machineNo;
        int processingTime;
};

class Operation{
    public:
        int jobId;
        int opJobSeqId;
        int routeChoicesNo;
        vector<MachineData> routeData;
};

class OpSched{
    public:
        int startTime;
        int endTime;
};

bool sortinrev(const pair<double, int> &a, const pair<double, int> &b){
    return (a.first>b.first);
}

bool compareByMakespan(const Particle &a, const Particle &b)
{
    return a.makespan < b.makespan;
}

double cal_magnitude(double arr[], int totalNoOps){
    double value = 0.0;
    for(int i=0; i<2*totalNoOps; i++){
        value+= arr[i]*arr[i];
    }
    return value;
}

int evaluate_makespan(vector<double> pos, vector<Operation> operationData, int numberOfJobs, int numMachines, double x_delta, int max_numOps){
    int totalNoOps=operationData.size();
    vector<int> routeVector;//will contain the values of the route chosen for each operation
    int priorityVector[totalNoOps];//will contain the sequence in which the operations are to be executed
    //CONSTRUCTING THE ROUTE VECTOR
    for(int i=0; i<totalNoOps; i++){
        int x = round((operationData[i].routeChoicesNo-1)*(pos[i]+x_delta)*(1/(2*x_delta)) + 1);
        routeVector.push_back(x-1);
    }
    //CONSTRUCTING THE SEQUENCE PRIORITY VECTOR
    vector<pair<double, int>> intermediateVec;
    for(int i=0; i<totalNoOps; i++){
        intermediateVec.push_back(make_pair(pos[i+totalNoOps], i));
    }
    sort(intermediateVec.begin(), intermediateVec.end(), sortinrev);
    int q[numberOfJobs]={0};
    for(int i=0; i<totalNoOps; i++){
        int job_Id=operationData[intermediateVec[i].second].jobId;
        int job_seq_Id=operationData[intermediateVec[i].second].opJobSeqId;
        q[job_Id]++;
        if(q[job_Id]-1 == job_seq_Id){
            priorityVector[i]=intermediateVec[i].second;
        }
        else{
            int placeCount=0;
            for(int j=0; j<totalNoOps; j++){
                if(operationData[intermediateVec[j].second].jobId==job_Id){
                    placeCount++;
                    if(placeCount-1==job_seq_Id){
                        priorityVector[j]=intermediateVec[i].second;
                    }
                }
            }
        }
    }

    //PRINTING ROUTE AND PRIORITY VECTOR
    /*
    cout<<"Route vector=";
    for(int i=0; i<totalNoOps; i++){
        cout<<routeVector[i]<<" ";
    }
    cout<<endl;
    cout<<"Priority Vector=";
    for(int i=0; i<totalNoOps; i++){
        cout<<priorityVector[i]<<" ";
    }
    cout<<endl;
    */

    //CALCULATING MAKESPAN USING ROUTE VECTOR AND SEQUENCE PRIORITY VECTOR
    int makespan=INT_MIN;
    vector<int> machineFreeTime(numMachines+1, 0);
    OpSched opMatrix[numberOfJobs][max_numOps];
    for(int i=0; i<totalNoOps; i++){
        int jobNo = operationData[priorityVector[i]].jobId;
        int opSeqId = operationData[priorityVector[i]].opJobSeqId;
        int routeChosen_i = routeVector[priorityVector[i]];  
        MachineData t = operationData[priorityVector[i]].routeData[routeChosen_i];
        if(opSeqId == 0){
            opMatrix[jobNo][opSeqId].startTime = machineFreeTime[t.machineNo];
            opMatrix[jobNo][opSeqId].endTime = opMatrix[jobNo][opSeqId].startTime + t.processingTime;
        }
        else{
            opMatrix[jobNo][opSeqId].startTime = max(machineFreeTime[t.machineNo], opMatrix[jobNo][opSeqId-1].endTime);
            opMatrix[jobNo][opSeqId].endTime = opMatrix[jobNo][opSeqId].startTime + t.processingTime;
        }
        machineFreeTime[t.machineNo]=opMatrix[jobNo][opSeqId].endTime;
        if(opMatrix[jobNo][opSeqId].endTime>makespan){
            makespan = opMatrix[jobNo][opSeqId].endTime;
        }
    }
    return makespan;
}

int main(){
    clock_t t;
    t=clock();
    ifstream inFile;
    inFile.open("kacem3.txt");
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    int numberOfJobs, totalNoOps, numMachines;
    inFile>>numberOfJobs;
    inFile>>numMachines;
    vector<Operation> operationData;
    int numOps;
    int max_numOps=INT_MIN;
    for(int i=0; i<numberOfJobs; i++){
        inFile>>numOps;
        max_numOps = max(max_numOps, numOps);
        for(int j=0; j<numOps; j++){
            Operation op1;
            op1.jobId=i;
            op1.opJobSeqId=j;
            inFile>>op1.routeChoicesNo;
            for(int k=0; k<op1.routeChoicesNo; k++){
                MachineData md1;
                inFile>>md1.machineNo;
                inFile>>md1.processingTime;
                op1.routeData.push_back(md1);
            }
            operationData.push_back(op1);
        }
    }
    inFile.close();

    //PRINTING THE INPUT VALUES FOR CHECKING
    /*
    for(int i=0; i<operationData.size(); i++){
        cout<<"Job ID="<<operationData[i].jobId<<" OperationSeqId="<<operationData[i].opJobSeqId<<endl;
        cout<<"RouteChoicesNo="<<operationData[i].routeChoicesNo<<endl;
        for(int j=0; j<operationData[i].routeData.size(); j++){
            cout<<"Machine No="<<operationData[i].routeData[j].machineNo<<" ProcessingTime="<<operationData[i].routeData[j].processingTime<<endl;;
        }
    }
    */

    /************************************************************************************************
        VARIABLES COMMON TO BOTH GREYWOLF AND WHALE
    *************************************************************************************************/
    totalNoOps=operationData.size();
    int populationSize = 50;
    double x_delta = 1.0;
    unsigned s = 2;
    default_random_engine generator (s);
    uniform_real_distribution<double> distribution(-1*x_delta,x_delta);
    uniform_real_distribution<double> unitDistribution(0,1);
    uniform_real_distribution<double> symmetricUnitDistribution(-1, 1);
    int num_iter=100;

    /*************************************************************************************************
        VARIABLES FOR USAGE WITH GREYWOLF OPTIMIZATION
    **************************************************************************************************/

    //GENERATION OF THE INITIAL RANDOM POPULATION
    vector<Particle> populationGrey;
    for(int i=0; i<populationSize; i++){
        Particle p1;
        for(int j=0; j<2*totalNoOps; j++){
            p1.pos.push_back(distribution(generator));
        }
        populationGrey.push_back(p1);
    }

    //INITIALIZE a, A, C
    double aGrey = 2.0;
    double AGrey[2*totalNoOps];
    double CGrey[2*totalNoOps];
    for(int i=0; i<2*totalNoOps; i++){
        AGrey[i] = 2*aGrey*(unitDistribution(generator)) - aGrey;
        CGrey[i] = 2*(unitDistribution(generator));
    }
    double x1[2*totalNoOps];
    double x2[2*totalNoOps];
    double x3[2*totalNoOps];


    /*************************************************************************************************
        VARIABLES FOR USAGE WITH WHALE OPTIMIZATION
    **************************************************************************************************/
    Particle p_best_W;
    double aWhale, l, p;
    double AWhale[2*totalNoOps];
    double CWhale[2*totalNoOps];
    //POPULATION INITIALIZATION for whale
    vector<Particle> populationWhale;
    for(int i=0; i<populationSize; i++){
        Particle p1;
        for(int j=0; j<2*totalNoOps; j++){
            p1.pos.push_back(distribution(generator));
        }
        populationWhale.push_back(p1);
    }

    //CALCULATE THE MAKESPAN FOR EACH PARTICLE in populationGrey
    for(int i=0; i<populationSize; i++){
        populationGrey[i].makespan = evaluate_makespan(populationGrey[i].pos, operationData, numberOfJobs, numMachines, x_delta, max_numOps);
    }

    //CALCULATE THE MAKESPAN FOR EACH PARTICLE in populationWhale
    for(int i=0; i<populationSize; i++){
        populationWhale[i].makespan = evaluate_makespan(populationWhale[i].pos, operationData, numberOfJobs, numMachines, x_delta, max_numOps);
    }
    
    //IDENTIFY THE THREE BEST SEARCH AGENTS-->for grey wolf
    sort(populationGrey.begin(), populationGrey.end(), compareByMakespan);
    Particle alpha;
    alpha.makespan = populationGrey[0].makespan;
    for(int i=0; i<2*totalNoOps; i++){
        alpha.pos.push_back(populationGrey[0].pos[i]);
    }
    Particle beta;
    beta.makespan = populationGrey[1].makespan;
    for(int i=0; i<2*totalNoOps; i++){
        beta.pos.push_back(populationGrey[1].pos[i]);
    }
    Particle delta;
    delta.makespan = populationGrey[2].makespan;
    for(int i=0; i<2*totalNoOps; i++){
        delta.pos.push_back(populationGrey[2].pos[i]);
    }

    //IDENTIFY THE BEST SEARCH AGENT-->for whale
    sort(populationWhale.begin(), populationWhale.end(), compareByMakespan);
    p_best_W.makespan = populationWhale[0].makespan;
    for(int i=0; i<2*totalNoOps; i++){
        p_best_W.pos.push_back(populationWhale[0].pos[i]);
    }

    //OPTIMIZATION PHASE 1 BEGINS
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                //SECTION CORRESPONDING TO GREY WOLF OPTIMIZATION
                for(int kg=1; kg<=num_iter/2; kg++){
                //UPDATING POSITION OF EACH SEARCH AGENT USING ALPHA, BETA, DELTA
                    for(int i=0; i<populationSize; i++){
                        for(int j=0; j<2*totalNoOps; j++){
                            x1[j] = alpha.pos[j] - AGrey[j]*(CGrey[j]*alpha.pos[j] - populationGrey[i].pos[j]);
                            x2[j] = beta.pos[j] - AGrey[j]*(CGrey[j]*beta.pos[j] - populationGrey[i].pos[j]);
                            x3[j] = delta.pos[j] - AGrey[j]*(CGrey[j]*delta.pos[j] - populationGrey[i].pos[j]);
                            populationGrey[i].pos[j] = (x1[j]+x2[j]+x3[j])/3;
                            if(populationGrey[i].pos[j]<-1*x_delta){
                                populationGrey[i].pos[j] = -1*x_delta;
                            }
                            if(populationGrey[i].pos[j]>x_delta){
                                populationGrey[i].pos[j] = x_delta;
                            }
                        }
                    }
                    //UPDATE a, A, C
                    aGrey = 2.0 - (2*kg)/num_iter;
                    for(int i=0; i<2*totalNoOps; i++){
                        AGrey[i] = 2*aGrey*(unitDistribution(generator)) - aGrey;
                        CGrey[i] = 2*(unitDistribution(generator));
                    }
                    //CALCULATE THE FITNESS OF ALL SEARCH AGENTS
                    for(int i=0; i<populationSize; i++){
                        populationGrey[i].makespan = evaluate_makespan(populationGrey[i].pos, operationData, numberOfJobs, numMachines, x_delta, max_numOps);
                    }
                    //UPDATE ALPHA, BETA, DELTA
                    sort(populationGrey.begin(), populationGrey.end(), compareByMakespan);
                    if(alpha.makespan > populationGrey[0].makespan){
                        alpha.makespan = populationGrey[0].makespan;
                        for(int i=0; i<2*totalNoOps; i++){
                            alpha.pos[i] = populationGrey[0].pos[i];
                        }
                    }
                    if(populationGrey[1].makespan < beta.makespan && populationGrey[1].makespan > alpha.makespan){
                        beta.makespan = populationGrey[1].makespan;
                        for(int i=0; i<2*totalNoOps; i++){
                            beta.pos[i] = populationGrey[1].pos[i];
                        }
                    }
                    if(populationGrey[2].makespan > alpha.makespan && populationGrey[2].makespan > beta.makespan && populationGrey[2].makespan < delta.makespan){
                        delta.makespan = populationGrey[2].makespan;
                        for(int i=0; i<2*totalNoOps; i++){
                            delta.pos[i] = populationGrey[2].pos[i];
                        }
                    }
                    //cout<<"G1 kg="<<kg<<" best_makespan="<<alpha.makespan<<endl;
                }
            }
            #pragma omp section
            {
                //SECTION CORRESPONDING TO WHALE OPTIMIZATION
                for(int kw=1; kw<=num_iter/2; kw++){
                    for(int i=0; i<populationSize; i++){
                        //UPDATE a, A, C, l and p
                        aWhale =2.0 -  2.0*kw/num_iter;
                        for(int j=0; j<2*totalNoOps; j++){
                            AWhale[j] = 2*aWhale*(unitDistribution(generator)) - aWhale;
                            CWhale[j] = 2*(unitDistribution(generator));
                        }
                        l = symmetricUnitDistribution(generator);
                        p = unitDistribution(generator);

                        if(p<0.5){
                            double mag_A = cal_magnitude(AWhale, totalNoOps);
                            if(mag_A<1){
                                for(int j=0; j<2*totalNoOps; j++){
                                    populationWhale[i].pos[j] = p_best_W.pos[j] - AWhale[j]*(CWhale[j]*p_best_W.pos[j] - populationWhale[i].pos[j]);
                                }
                            }
                            else{
                                int r_i = rand()%populationSize;
                                for(int j=0; j<2*totalNoOps; j++){
                                    populationWhale[i].pos[j] = populationWhale[r_i].pos[j] - AWhale[j]*(CWhale[j]*populationWhale[r_i].pos[j] - populationWhale[i].pos[j]);
                                }
                            }
                        }
                        else{
                            for(int j=0; j<2*totalNoOps; j++){
                                populationWhale[i].pos[j] = (p_best_W.pos[j] - populationWhale[i].pos[j])*exp(0.206*l)*cos(2*3.141592*l) + p_best_W.pos[j];
                            }
                        }

                        //CHECKING FOR EXCEED OF BOUNDS
                        for(int j=0; j<2*totalNoOps; j++){
                            if(populationWhale[i].pos[j] < -1*x_delta){
                                populationWhale[i].pos[j] = -1*x_delta;
                            }
                            if(populationWhale[i].pos[j] > x_delta){
                                populationWhale[i].pos[j] = x_delta;
                            }
                        }
                        populationWhale[i].makespan = evaluate_makespan(populationWhale[i].pos, operationData, numberOfJobs, numMachines, x_delta, max_numOps);
                    }
                    sort(populationWhale.begin(), populationWhale.end(), compareByMakespan);
                    if(populationWhale[0].makespan<p_best_W.makespan){
                        p_best_W.makespan = populationWhale[0].makespan;
                        for(int j=0; j<2*totalNoOps; j++){
                            p_best_W.pos[j] = populationWhale[0].pos[j];
                        }
                    }
                    //cout<<"W1 kw="<<kw<<" best_makespan="<<p_best_W.makespan<<endl;
                }
            }
        }
        #pragma omp barrier
    }

    //IDENTIFY THE THREE BEST SEARCH AGENTS-->for grey wolf
    alpha.makespan = populationWhale[0].makespan;
    for(int i=0; i<2*totalNoOps; i++){
        alpha.pos[i] = populationWhale[0].pos[i];
    }
    beta.makespan = populationWhale[1].makespan;
    for(int i=0; i<2*totalNoOps; i++){
        beta.pos[i] = populationWhale[1].pos[i];
    }
    delta.makespan = populationWhale[2].makespan;
    for(int i=0; i<2*totalNoOps; i++){
        delta.pos[i] = populationWhale[2].pos[i];
    }

    //IDENTIFY THE BEST SEARCH AGENT-->for whale
    p_best_W.makespan = populationGrey[0].makespan;
    for(int i=0; i<2*totalNoOps; i++){
        p_best_W.pos[i] = populationGrey[0].pos[i];
    }

    //OPTIMIZATION PHASE 2 BEGINS
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                //SECTION CORRESPONDING TO GREY WOLF OPTIMIZATION
                for(int kg=1; kg<=num_iter/2; kg++){
                //UPDATING POSITION OF EACH SEARCH AGENT USING ALPHA, BETA, DELTA
                    for(int i=0; i<populationSize; i++){
                        for(int j=0; j<2*totalNoOps; j++){
                            x1[j] = alpha.pos[j] - AGrey[j]*(CGrey[j]*alpha.pos[j] - populationWhale[i].pos[j]);
                            x2[j] = beta.pos[j] - AGrey[j]*(CGrey[j]*beta.pos[j] - populationWhale[i].pos[j]);
                            x3[j] = delta.pos[j] - AGrey[j]*(CGrey[j]*delta.pos[j] - populationWhale[i].pos[j]);
                            populationWhale[i].pos[j] = (x1[j]+x2[j]+x3[j])/3;
                            if(populationWhale[i].pos[j]<-1*x_delta){
                                populationWhale[i].pos[j] = -1*x_delta;
                            }
                            if(populationWhale[i].pos[j]>x_delta){
                                populationWhale[i].pos[j] = x_delta;
                            }
                        }
                    }
                    //UPDATE a, A, C
                    aGrey = 2.0 - (2*kg)/num_iter;
                    for(int i=0; i<2*totalNoOps; i++){
                        AGrey[i] = 2*aGrey*(unitDistribution(generator)) - aGrey;
                        CGrey[i] = 2*(unitDistribution(generator));
                    }
                    //CALCULATE THE FITNESS OF ALL SEARCH AGENTS
                    for(int i=0; i<populationSize; i++){
                        populationWhale[i].makespan = evaluate_makespan(populationWhale[i].pos, operationData, numberOfJobs, numMachines, x_delta, max_numOps);
                    }
                    //UPDATE ALPHA, BETA, DELTA
                    sort(populationWhale.begin(), populationWhale.end(), compareByMakespan);
                    if(alpha.makespan > populationWhale[0].makespan){
                        alpha.makespan = populationWhale[0].makespan;
                        for(int i=0; i<2*totalNoOps; i++){
                            alpha.pos[i] = populationWhale[0].pos[i];
                        }
                    }
                    if(populationWhale[1].makespan < beta.makespan && populationWhale[1].makespan > alpha.makespan){
                        beta.makespan = populationWhale[1].makespan;
                        for(int i=0; i<2*totalNoOps; i++){
                            beta.pos[i] = populationWhale[1].pos[i];
                        }
                    }
                    if(populationWhale[2].makespan > alpha.makespan && populationWhale[2].makespan > beta.makespan && populationWhale[2].makespan < delta.makespan){
                        delta.makespan = populationWhale[2].makespan;
                        for(int i=0; i<2*totalNoOps; i++){
                            delta.pos[i] = populationWhale[2].pos[i];
                        }
                    }
                    //cout<<"G2 kg="<<kg<<" best_makespan="<<alpha.makespan<<endl;
                }
            }
            #pragma omp section
            {
                //SECTION CORRESPONDING TO WHALE OPTIMIZATION
                for(int kw=1; kw<=num_iter/2; kw++){
                    for(int i=0; i<populationSize; i++){
                        //UPDATE a, A, C, l and p
                        aWhale =2.0 -  2.0*kw/num_iter;
                        for(int j=0; j<2*totalNoOps; j++){
                            AWhale[j] = 2*aWhale*(unitDistribution(generator)) - aWhale;
                            CWhale[j] = 2*(unitDistribution(generator));
                        }
                        l = symmetricUnitDistribution(generator);
                        p = unitDistribution(generator);

                        if(p<0.5){
                            double mag_A = cal_magnitude(AWhale, totalNoOps);
                            if(mag_A<1){
                                for(int j=0; j<2*totalNoOps; j++){
                                    populationGrey[i].pos[j] = p_best_W.pos[j] - AWhale[j]*(CWhale[j]*p_best_W.pos[j] - populationGrey[i].pos[j]);
                                }
                            }
                            else{
                                int r_i = rand()%populationSize;
                                for(int j=0; j<2*totalNoOps; j++){
                                    populationGrey[i].pos[j] = populationGrey[r_i].pos[j] - AWhale[j]*(CWhale[j]*populationGrey[r_i].pos[j] - populationGrey[i].pos[j]);
                                }
                            }
                        }
                        else{
                            for(int j=0; j<2*totalNoOps; j++){
                                populationGrey[i].pos[j] = (p_best_W.pos[j] - populationGrey[i].pos[j])*exp(0.206*l)*cos(2*3.141592*l) + p_best_W.pos[j];
                            }
                        }

                        //CHECKING FOR EXCEED OF BOUNDS
                        for(int j=0; j<2*totalNoOps; j++){
                            if(populationGrey[i].pos[j] < -1*x_delta){
                                populationGrey[i].pos[j] = -1*x_delta;
                            }
                            if(populationGrey[i].pos[j] > x_delta){
                                populationGrey[i].pos[j] = x_delta;
                            }
                        }
                        populationGrey[i].makespan = evaluate_makespan(populationGrey[i].pos, operationData, numberOfJobs, numMachines, x_delta, max_numOps);
                    }
                    sort(populationGrey.begin(), populationGrey.end(), compareByMakespan);
                    if(populationGrey[0].makespan<p_best_W.makespan){
                        p_best_W.makespan = populationGrey[0].makespan;
                        for(int j=0; j<2*totalNoOps; j++){
                            p_best_W.pos[j] = populationGrey[0].pos[j];
                        }
                    }
                    //cout<<"W2 kw="<<kw<<" best_makespan="<<p_best_W.makespan<<endl;
                }
            }
        }
        #pragma omp barrier
    }

    //COMPARISON OF THE RESULTS OBTAINED FROM THE TWO OPTIMIZATION PARADIGMS
    if(alpha.makespan < p_best_W.makespan){
        cout<<"Best Solution:"<<"\n";
        cout<<"Makespan:"<<alpha.makespan<<"\n";
        cout<<endl;
        /*

        vector<int> routeVector;//will contain the values of the route chosen for each operation
        int priorityVector[totalNoOps];
        //CONSTRUCTING THE ROUTE VECTOR
        for(int i=0; i<totalNoOps; i++){
            int x = round((operationData[i].routeChoicesNo-1)*(alpha.pos[i]+x_delta)*(1/(2*x_delta)) + 1);
            routeVector.push_back(x-1);
        }
        //CONSTRUCTING THE SEQUENCE PRIORITY VECTOR
        vector<pair<double, int>> intermediateVec;
        for(int i=0; i<totalNoOps; i++){
            intermediateVec.push_back(make_pair(alpha.pos[i+totalNoOps], i));
        }
        sort(intermediateVec.begin(), intermediateVec.end(), sortinrev);
        int q[numberOfJobs]={0};
        for(int i=0; i<totalNoOps; i++){
            int job_Id=operationData[intermediateVec[i].second].jobId;
            int job_seq_Id=operationData[intermediateVec[i].second].opJobSeqId;
            q[job_Id]++;
            if(q[job_Id]-1 == job_seq_Id){
                priorityVector[i]=intermediateVec[i].second;
            }
            else{
                int placeCount=0;
                for(int j=0; j<totalNoOps; j++){
                    if(operationData[intermediateVec[j].second].jobId==job_Id){
                        placeCount++;
                        if(placeCount-1==job_seq_Id){
                            priorityVector[j]=intermediateVec[i].second;
                        }
                    }
                }
            }
        }
            //PRINTING ROUTE AND PRIORITY VECTOR
        
        cout<<"Route vector=";
        for(int i=0; i<totalNoOps; i++){
            cout<<routeVector[i]<<" ";
        }
        cout<<endl;
        cout<<"Priority Vector=";
        for(int i=0; i<totalNoOps; i++){
            cout<<priorityVector[i]<<" ";
        }
        cout<<endl;
        */
        


    }
    else{
        cout<<"Best Solution:"<<"\n";
        cout<<"Makespan:"<<p_best_W.makespan<<"\n";
        cout<<endl;

        /*
        vector<int> routeVector;//will contain the values of the route chosen for each operation
        int priorityVector[totalNoOps];
        //CONSTRUCTING THE ROUTE VECTOR
        for(int i=0; i<totalNoOps; i++){
            int x = round((operationData[i].routeChoicesNo-1)*(p_best_W.pos[i]+x_delta)*(1/(2*x_delta)) + 1);
            routeVector.push_back(x-1);
        }
        //CONSTRUCTING THE SEQUENCE PRIORITY VECTOR
        vector<pair<double, int>> intermediateVec;
        for(int i=0; i<totalNoOps; i++){
            intermediateVec.push_back(make_pair(p_best_W.pos[i+totalNoOps], i));
        }
        sort(intermediateVec.begin(), intermediateVec.end(), sortinrev);
        int q[numberOfJobs]={0};
        for(int i=0; i<totalNoOps; i++){
            int job_Id=operationData[intermediateVec[i].second].jobId;
            int job_seq_Id=operationData[intermediateVec[i].second].opJobSeqId;
            q[job_Id]++;
            if(q[job_Id]-1 == job_seq_Id){
                priorityVector[i]=intermediateVec[i].second;
            }
            else{
                int placeCount=0;
                for(int j=0; j<totalNoOps; j++){
                    if(operationData[intermediateVec[j].second].jobId==job_Id){
                        placeCount++;
                        if(placeCount-1==job_seq_Id){
                            priorityVector[j]=intermediateVec[i].second;
                        }
                    }
                }
            }
        }
        t = clock()-t;
        double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
        cout<<"time taken ="<<time_taken<<endl;

        //PRINTING ROUTE AND PRIORITY VECTOR
        cout<<"Route vector=";
        for(int i=0; i<totalNoOps; i++){
            cout<<routeVector[i]<<" ";
        }
        cout<<endl;
        cout<<"Priority Vector=";
        for(int i=0; i<totalNoOps; i++){
            cout<<priorityVector[i]<<" ";
        }
        cout<<endl;
        */
        

    }
    t = clock()-t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    cout<<"time taken ="<<time_taken<<endl;

    return 0;
}