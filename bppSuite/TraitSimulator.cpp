// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// From bpp-core:
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h> 
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceTools.h>

// From bpp-phyl
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Simulation/DetailedSiteSimulator.h>
#include <Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>

using namespace bpp;
using namespace std;

typedef vector<vector<double>> VVDouble;
typedef vector<double> VDouble;
typedef unsigned int uint;

#define STATE "state"



void giveNamesToInternalNodes(Tree* tree)
{
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);
    vector<Node*> nodes = ttree->getNodes();
    for (size_t i=0; i<nodes.size(); ++i) {
        if (!nodes[i]->hasName())
            nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
    }  
}

Tree* processTree(BppApplication* bppml)
{
    Tree* tree = PhylogeneticsApplicationTools::getTree(bppml->getParams());
    giveNamesToInternalNodes(tree);
    TreeTemplate<Node> ttree(*tree);
    vector<Node *> nodes = ttree.getNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      if(nodes[i]->isLeaf())
        nodes[i]->setName("Leaf" + TextTools::toString(nodes[i]->getId()));
      else
        nodes[i]->setName("Internal" + TextTools::toString(nodes[i]->getId()));
      ApplicationTools::displayResult(nodes[i]->getName(), TextTools::toString(nodes[i]->getId())); 
    }
    
    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml->getParams(), "Input", "", true, 1); // process tree branches lengths
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs); // this line process cmdName - which dictates weather the tree should be processed as is, or ultrameterized
    if (cmdName == "Clock") // if needed, turn the tree into an ultrametric tree
    {
      TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
    }
    return tree;
}

int getNodeState(const Node* node)
{
    return (dynamic_cast<const BppInteger*>(node->getNodeProperty(STATE)))->getValue(); // exception on root on the true history - why didn't the root recieve a state?
}

void setNodeState(Node* node, size_t state)
{
    BppInteger* stateProperty = new BppInteger(static_cast<int>(state));
    node->setNodeProperty(STATE, *stateProperty);
    delete stateProperty;
}

// since the function is fitted to true history, the mutation path may exceed the original branch length, in which case the dutration until the last event should be reduced to fit the original branch length
void updateBranchMapping(Node* son, const MutationPath& branchMapping, size_t initial_count)
{ 
    static size_t nodesCounter_ = initial_count;
    const vector<size_t> states = branchMapping.getStates();
    const VDouble times = branchMapping.getTimes();
    Node* curNode = son; // builds the branch history bottom to top (from last event to first event)
    double origBranchLength = son->getDistanceToFather();
    Node* nextNode;
    int eventsNum = static_cast<int>(branchMapping.getNumberOfEvents());
    if (eventsNum == 0)  // if there are no events -> return nothing
    {
        return;
    }
    else
    {
        // update the branch length of the node to be as per the time of the last event
        double duration = 0;
        for (int i=eventsNum-1; i>-1; --i) 
        { // add a new node to represent the transition
            nodesCounter_ += 1;
            const string name = "_mappingInternal" + TextTools::toString(nodesCounter_) + "_";
            nextNode = new Node(static_cast<int>(nodesCounter_), name); // nodes counter isn't incremented properly
            if (states[i] == 0)
            {
                setNodeState(nextNode,1);
            }
            else 
            {
                setNodeState(nextNode,0);
            }
            if (i == 0)
            {
                duration = times[i];    
            }
            else
            {
                duration = times[i]-times[i-1]; // the duration is the gap between the two adjacent times
            }
            nextNode->setDistanceToFather(duration); // according to the simulator, the documented time is the accumulated time until the transition and NOT the time since the last transition
            
            // set the father to no longer be the father of curNode
            if (curNode->hasFather()) // as long as we haven't reached the root -> upate the father of the current node
            {
                Node* originalFather = curNode->getFather();
                originalFather->removeSon(curNode); // also removes originalFather as the father of curNode
                // set nextNode to be the new father of curNode
                curNode->setFather(nextNode); // also adds curNode to the sons of nextNode
                // set curNode's original father ot be the father of nextNode
                nextNode->setFather(originalFather); // also adds nextNode to the sons of originalFather - or not? make sure this is the father at all times
                curNode = nextNode;
            }
        }
        // if the sum of distances doesn't add up to the length of the original branch, add an event or change the duration of the last event, depending on the last assigned state
        double lastEventDuration = origBranchLength - times[eventsNum-1]; 
        if (lastEventDuration > 0)
        {
            son->setDistanceToFather(lastEventDuration);
        }
    }
    return;
}

void updateStatesInNodesNames(Tree* mapping)
{
    string label = "state";
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(mapping))->getNodes(); 
    for (int i=0; i < static_cast<int>(nodes.size()); i++) 
    {
        string name = nodes[i]->getName();
        int state = getNodeState(nodes[i]);
        nodes[i]->setName(name + "{" + TextTools::toString(state) + "}");
    }
}

// translate the simulation data into a tree format, similar to the one given by the stochastic mappings
// in the results, there is a mutations paths field that can be returned for each node getMutationPath(int nodeId)
// then, it is enough to use the function I built in StochasticMapping
Tree* extractTrueHistory(RASiteSimulationResult* simulationData, Tree* baseTree)
{
    Tree* history = baseTree->clone();
    giveNamesToInternalNodes(history);
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(history))->getNodes();
    vector<Node*> debugNodes;
    size_t state;
    // set the ancestral states for all the nodes
    for (size_t j=0; j<nodes.size(); ++j)
    {
        // bypass bug in DetailedSiteSimulator.h line 127 - for root ancestral state, need to access ancestralStates[0] but can't because its index in indexes_ is 0 -> so line 127 end up accessing index 1
        // to bypass, call function in line 125 that accesses the argument index directly
        string nodeName = nodes[j]->getName();
        if (nodes[j]->getId() == history->getRootId())
        {
            state = simulationData->getRootAncestralState();
        }
        else
        {
            state = simulationData->getAncestralState(nodes[j]->getId());
        }
        setNodeState(nodes[j], state);
    }
    for (size_t i=0; i<nodes.size(); ++i)
    {
        string name = nodes[i]->getName();
        if (nodes[i]->hasFather())
        {
            // extract the path that ends in the visited node and update it on the tree
            const MutationPath branchHistory = simulationData->getMutationPath(nodes[i]->getId());
            updateBranchMapping(nodes[i], branchHistory, nodes.size()-1);
        } 
    }
    return history;
}


void checkTrueHistory(Tree* trueHistory, Tree* baseTree)
{
    // make sure both trees have the same size
    VDouble baseTreeBranchLengths = dynamic_cast<TreeTemplate<Node>*>(baseTree)->getBranchLengths();
    double baseTreeSize = 0;
    for (size_t b=0; b<baseTreeBranchLengths.size(); b++)
    {
        baseTreeSize = baseTreeSize + baseTreeBranchLengths[b];
    }

    VDouble trueHistoryBranchLengths = dynamic_cast<TreeTemplate<Node>*>(trueHistory)->getBranchLengths();
    double trueHistorySize = 0;
    for (size_t b=0; b < trueHistoryBranchLengths.size(); b++)
    {
        trueHistorySize = trueHistorySize + trueHistoryBranchLengths[b];
    }

    if (abs(baseTreeSize - trueHistorySize) > 0.01)
    {
        throw Exception("true history has different tbl from base tree");
    }

    string nodeName;
    string mappingNodeName = "mapping";
    vector<Node*> trueHistoryNodes = (dynamic_cast<TreeTemplate<Node>*>(trueHistory))->getNodes();
    
    // make sure that each son of a base node has the same state as the base node
    Node* son;
    int nodeState, sonState;
    for (size_t n=0; n<trueHistoryNodes.size(); n++)
    {
        nodeName = trueHistoryNodes[n]->getName();
        // check if the node is a base node (either baseInternal of leaf)
        if (nodeName.find(mappingNodeName) == std::string::npos)
        {
            nodeState = getNodeState(trueHistoryNodes[n]);
            size_t numberOfSons = trueHistoryNodes[n]->getNumberOfSons();
            for (size_t s=0; s<numberOfSons; s++)
            {
                son = trueHistoryNodes[n]->getSon(s);
                sonState = getNodeState(son);
                if (sonState != nodeState)
                {
                    // if the son is a mapping node, its state should be the same as the state of its father base node
                    if (son->getName().find(mappingNodeName) != std::string::npos)
                    {
                        throw Exception("son state doesn't match parent state at node: " + son->getName());
                    }
                }
            }
        }  
    }

    // make sure that each parent of a base node has the same state as the base node
    int fatherState;
    Node* father;
    for(size_t i=0; i<trueHistoryNodes.size(); i++)
    {
        nodeName = trueHistoryNodes[i]->getName();
        // check if the node is a base node (either baseInternal of leaf)
        if (trueHistoryNodes[i]->getId() != trueHistory->getRootId())
        {
            if (!(nodeName.find(mappingNodeName) != std::string::npos))
            {
                nodeState = getNodeState(trueHistoryNodes[i]);
                father = trueHistoryNodes[i]-> getFather();
                fatherState = getNodeState(father);
                if (fatherState == nodeState)
                {
                    // if the parent is a mappimg node, its state must be differnt from the state of its base tree son
                    if (father->getName().find(mappingNodeName) != std::string::npos)
                    {
                        throw Exception("father state matches parent state at node: " + father->getName());
                    }
                    
                }
            }
        }  
    }
}

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char** argv)
{
  cout << "************************************************************************************************" << endl;
  cout << "*                                     Trait simulator                                          *" << endl;
  cout << "************************************************************************************************" << endl;
  cout << endl;

    //long mySeed = 5698865; //RandomTools::giveRandomNumberBetweenZeroAndEntry(10000000);
    //cout << "mySeed : " << mySeed << endl;
    //RandomTools::setSeed(mySeed);

  try
  {
    /* process input from params file */
    BppApplication simulationParams(args, argv, "simulationParams");

    // extract from the parameters files the simulation parameters
    double mu = ApplicationTools::getDoubleParameter("character_model.mu", simulationParams.getParams(), 10., "", true, 2);
    double pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", simulationParams.getParams(), 0.5, "", true, 2);;

    // process the base tree
    Tree* tree = processTree(&simulationParams);

    // create a binary model
    const BinaryAlphabet* alphabet = new BinaryAlphabet();
    TwoParameterBinarySubstitutionModel* charModel = new TwoParameterBinarySubstitutionModel(alphabet,mu, pi0); // second arguent stands for mu
    DiscreteDistribution* rdist = new ConstantRateDistribution();
    vector<string> seqNames = tree->getLeavesNames();
    VectorSiteContainer charData(seqNames, alphabet);
     
    // simulate character history using a simulator over a simple binary model
    NonHomogeneousSequenceSimulator* charSimulator = new NonHomogeneousSequenceSimulator(charModel, rdist, tree);
    
    
    RASiteSimulationResult* charResult = charSimulator->dSimulateSite(); // dSimulateSite is a faulty function! but I need the elavorated result to simulate!
    unique_ptr<Site> charSite(charResult->getSite(*charSimulator->getSubstitutionModelSet()->getModel(0)));
    charSite->setPosition(0);
    charData.addSite(*charSite, false);
    Tree* trueHistory = extractTrueHistory(charResult, tree);

    // check the legality of the true history
    checkTrueHistory(trueHistory, tree);

    // write the simulate character data
    std::map<std::string, std::string> writeParams;
    writeParams["output.sequence.file"] = ApplicationTools::getAFilePath("character_history.seq_path", simulationParams.getParams(), false, false, "", true, "none", 1);
    writeParams["output.sequence.format"] = "Fasta";
    SequenceApplicationTools::writeSequenceFile(*(dynamic_cast<SequenceContainer*>(&charData)), writeParams, "", false, 1);

    // write the true history to a file before deleting it
    writeParams["output.tree.file"]  = ApplicationTools::getAFilePath("character_history.tree_path", simulationParams.getParams(), false, false, "", true, "none", 1);
    writeParams["output.tree.format"] = "Newick";
    // now, add the label of each internal node in the updated history to the node's name
    updateStatesInNodesNames(trueHistory);
    PhylogeneticsApplicationTools::writeTree(*trueHistory, writeParams);
    delete trueHistory;
    
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}