#include <bits/stdc++.h>
#include <random>
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
    ifstream inFile;
    inFile.open("datafile.txt");
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

    //GENERATION OF THE INITIAL RANDOM POPULATION
    totalNoOps=operationData.size();
    int populationSize = 10;
    vector<Particle> population;
    double x_delta = 3.0;
    unsigned s = 2;
    default_random_engine generator (s);
    uniform_real_distribution<double> distribution(-1*x_delta,x_delta);
    for(int i=0; i<populationSize; i++){
        Particle p1;
        for(int j=0; j<2*totalNoOps; j++){
            p1.pos.push_back(distribution(generator));
        }
        population.push_back(p1);
    }

    //PRINTING THE POPULATION VECTORS
    /*
    for(int i=0; i<populationSize; i++){
        for(int j=0; j<2*totalNoOps; j++){
            cout<<population[i].pos[j]<<" ";
        }
        cout<<endl;
    }
    */

    //CALCULATE THE MAKESPAN FOR EACH PARTICLE
    for(int i=0; i<populationSize; i++){
        population[i].makespan = evaluate_makespan(population[i].pos, operationData, numberOfJobs, numMachines, x_delta, max_numOps);
        //cout<<"makespan at i="<<i<<" is "<<population[i].makespan<<endl;
    }
    

    //FIND OUT ALPHA BETA DELTA FOR THE CURRENT POPULATION


    //OPTIMIZATION LOOP BEGINS

    //UPDATE THE POSITIONS OF ALL PARTICLES IN THE POPULATION USING ALPHA BETA AND DELTA

    //CHECK FOR TERMINATION CONDITION AND REPEAT

    //PRINT ALPHA



    return 0;
}