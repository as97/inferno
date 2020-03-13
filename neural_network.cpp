#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

using namespace std;

struct connection
{
    double weight;
    double deltaWeight;
};


class Neuron;

typedef vector<Neuron> Layer;

//******************Class Neuron*********************
class Neuron
{
public:
    Neuron(unsigned numOutput, unsigned myIndex);
    void setoutVal(double val){ n_outVal = val;}
    double getoutVal(void) const { return n_outVal; }
    void feedforward(const Layer &prevLayer);

private:
    static double activationFunction(double x);
    static double activationDerivative(double x);
    static double randomWeight(void) { return rand()/double(RAND_MAX); }
    double n_outVal;
    unsigned n_myIndex;
    vector<connection> n_outWeight;
};
double Neuron::activationFunction(double x){//Hyperbolic func
    return tanh(x);
}
double Neuron::activationDerivative(double x){
    return 1 - ( tanh(x))*tanh(x);
}
void Neuron::feedforward(const Layer &prevLayer){
    double sum = 0.0;
    for (unsigned n = 0; n < prevLayer.size(); n++){
         sum = sum + prevLayer[n].getoutVal()*prevLayer[n].n_outWeight[n_myIndex].weight;
    }
    n_outVal = Neuron::activationFunction(sum);
}
Neuron::Neuron(unsigned numOutput, unsigned myIndex){
    for (unsigned c = 0; c < numOutput; c++){
        n_outWeight.push_back(connection());
        n_outWeight.back().weight = randomWeight();
        n_myIndex= myIndex;
    }
}


//******************Class Net************************
class Net
{
public:
    Net(const vector<unsigned> &topology);//constructor
    void feedForward(const vector<double> &inputVal);//empty body declaration of member function
    void backProp(const vector<double> &targetVal);
    void getResult(vector<double> &resultVal) const {};//const correctness: member function that reads output values and put them in the container and doesn't modify the object of this class at all

private:
    vector<Layer> n_layer; //m_layer[layerNum][neuronNum]
    double n_error;
    double n_recentAvgError;
    double smoothingFactor;
};
void Net::backProp(const vector<double> &targetVal)
{
    // Calculate overall net error (RMS of output neuron error)
    Layer &outputLayer = n_layer.back();
    n_error = 0.0;
    for(unsigned n = 0; n < outputLayer.size(); n++){
        double delta = targetVal[n] - outputLayer[n].getoutVal();
        n_error = delta*delta;
    }
    n_error /= outputLayer.size() - 1;
    n_error = sqrt(n_error);

    // Implement a recent average measurement
    n_recentAvgError = (n_recentAvgError*smoothingFactor + n_error)/(smoothingFactor+1);

    // Calculate output layer gradient

    // Calculate gradient on hidden layer

    // For all  layer from outputs to first hidden layer

    // Update connection weights

}
void Net::feedForward(const vector<double> &inputVal){
    cout << "Input size = " << inputVal.size() << " Input layer size = " << n_layer[0].size() <<endl;
    assert(inputVal.size() == n_layer[0].size());

    //assign inputval to input neurons
    for( unsigned i = 0; i < inputVal.size(); i++ ){
        n_layer[0][i].setoutVal(inputVal[i]);
    }

    //Forward propagation
    for(unsigned numLayer = 1; numLayer < n_layer.size(); numLayer++ ){
        Layer &prevLayer = n_layer[numLayer-1]; 
        for(unsigned n = 0; n < n_layer[numLayer].size(); n++){
            n_layer[numLayer][n].feedforward(prevLayer);
        }
    }
}

Net::Net(const vector<unsigned> &topology)
{
    unsigned numLayer = topology.size();
    for ( unsigned layerNum = 0; layerNum < numLayer; layerNum++ ){
        cout << "Building Layer "<< layerNum+1 << endl; 
        n_layer.push_back(Layer());//appends an object into the container(m_layer) and the value is the Layer constructor
        unsigned numOutput = layerNum == numLayer - 1 ? 0 : topology[layerNum + 1];

        for ( unsigned neuronNum = 0; neuronNum < topology[layerNum]; neuronNum++){
            n_layer.back().push_back(Neuron(numOutput, neuronNum));//.back() access the most recently appended object into the container(m_layer) and the value is the Neuron constructor
            cout << "\t Making Neuron " << neuronNum+1 << endl;
        }
    }
}

int main()
{
    //eg, (3,2,1)
    vector<unsigned> topology;//defines the topology of Neural net
    topology.push_back(3);
    topology.push_back(2);
    topology.push_back(1);

    Net myNet(topology);

    vector<double> inputVal;//variable dim arrray or vector of type double
    myNet.feedForward(inputVal);
    
    vector<double> targetVal;
    myNet.backProp(targetVal);

    vector<double> resultVal;
    myNet.getResult(resultVal);
}
