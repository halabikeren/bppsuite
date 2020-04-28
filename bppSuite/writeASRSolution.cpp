// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/BppOSequenceReaderFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentReaderFormat.h>
#include <Bpp/Seq/Io/BppOSequenceWriterFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>


// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequencySet/MvaFrequencySet.h>
#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOFrequencySetFormat.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/Likelihood/JointLikelihoodFunction.h>
//#include <Bpp/Phyl/AncestralStateReconstruction.h>
#include <Bpp/Phyl/Likelihood/MarginalAncestralStateReconstruction.h>

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <vector>
//#include <boost/math/distributions/chi_squared.hpp>


using namespace bpp;
using namespace std;

typedef vector<vector<double>> VVDouble;
typedef vector<double> VDouble;
typedef unsigned int uint;

#define STATE "state"

/******************************************************************************/
/**************************** Auxiliary functions *****************************/
/******************************************************************************/

const CodonAlphabet* getCodonAlphabet()
{
  map<string, string> alphabetSettings;
  alphabetSettings["alphabet"] = "Codon(letter=DNA)";
  const Alphabet* alphabet = SequenceApplicationTools::getAlphabet(alphabetSettings, "", false); 
  const CodonAlphabet* codonAlphabet = dynamic_cast<const CodonAlphabet*>(alphabet);
  return codonAlphabet;
}

/******************************************************************************/

VectorSiteContainer* processCharacterData(BppApplication* bppml, const BinaryAlphabet* alphabet)
{
    string charDataFilePath = ApplicationTools::getAFilePath("input.character.file", bppml->getParams(), true, true, "", true, "none", 1);
    string sequenceFormat = ApplicationTools::getStringParameter("input.character.format", bppml->getParams(), "Fasta()", "", true, 1);
    BppOAlignmentReaderFormat bppoReader(1);
    unique_ptr<IAlignment> iAln(bppoReader.read(sequenceFormat));
    map<string, string> args(bppoReader.getUnparsedArguments());
    ApplicationTools::displayResult("character data file ", charDataFilePath);
    ApplicationTools::displayResult("chatacter data format ", iAln->getFormatName());
    const SequenceContainer* charCont = iAln->readAlignment(charDataFilePath, alphabet);
    VectorSiteContainer* sites = new VectorSiteContainer(*dynamic_cast<const OrderedSequenceContainer*>(charCont));
    delete charCont;
    return sites;
}

/******************************************************************************/

VectorSiteContainer* processCodonAlignment(BppApplication* bppml, const CodonAlphabet* codonAlphabet)
{
    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(codonAlphabet, bppml->getParams()); // here, gaps will be converted to unknown characters
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml->getParams(), "", true, false); // convert the alignemnt to condensed format of unique sites
    delete allSites; // delete the non-condenced intance of the sequence data
    return sites;
}

/******************************************************************************/

void giveNamesToInternalNodes(Tree* tree)
{
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);
    vector<Node*> nodes = ttree->getNodes();
    for (size_t i=0; i<nodes.size(); ++i) {
        if (!nodes[i]->hasName())
            nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
    }  
}

/******************************************************************************/

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

/******************************************************************************/

TransitionModel* setCharacterModel(BppApplication* bppml, VectorSiteContainer* charData, const BinaryAlphabet* alphabet, Tree* tree)
{
    // create the model
    SubstitutionModel* model = new TwoParameterBinarySubstitutionModel(alphabet);
    return dynamic_cast<TransitionModel*>(model);
}

/******************************************************************************/

void setPartition(BppApplication* bppml, Tree* solution)
{

  // set assignment to branches
  vector <const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>*>(solution))->getNodes();
  string character0NodesIds, character1NodesIds = "";
  for (size_t i=0; i<nodes.size(); ++i)
  {
    int nodeState = dynamic_cast<const BppInteger*>(nodes[i]->getNodeProperty("state"))->getValue();
    if (nodeState == 0)
    {
      character0NodesIds = character0NodesIds + TextTools::toString(nodes[i]->getId()) + ","; 
    } 
    else 
    {
      character1NodesIds = character1NodesIds + TextTools::toString(nodes[i]->getId()) + ","; 
    }
  }

  string outputPath = ApplicationTools::getStringParameter("output.mp.partition", bppml->getParams(), "", "", true, 1);
  ofstream myfile (outputPath);
  if (myfile.is_open())
  {
    myfile << character0NodesIds.substr(0,character0NodesIds.length()-1) << "\n";
    myfile << character1NodesIds.substr(0,character1NodesIds.length()-1) << "\n";
  }
  myfile.close();
}

/******************************************************************************/

int getNodeState(const Node* node)
{
    return (dynamic_cast<const BppInteger*>(node->getNodeProperty("state")))->getValue();
}

/******************************************************************************/

void setNodeState(Node* node, size_t state)
{
    BppInteger* stateProperty = new BppInteger(static_cast<int>(state));
    node->setNodeProperty(STATE, *stateProperty);
    delete stateProperty;
}
/******************************************************************************/

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

/******************************************************************************/

Tree* computeASRMLSolution(DRHomogeneousTreeLikelihood* tl)
{
	const Tree& tree = tl->getTree();
	Tree* solution = tree.clone();
	vector<const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>&>(tree)).getNodes();
	AncestralStateReconstruction* asr = new MarginalAncestralStateReconstruction(tl);
	vector<size_t> nodeStates;
	for (size_t n=0; n<nodes.size(); ++n)
	{
		Node* copyNode = (dynamic_cast<TreeTemplate<Node>*>(solution))->getNode(nodes[n]->getId());
		if (nodes[n]->getId() != tree.getRootId())
		{
			nodeStates = asr->getAncestralStatesForNode(nodes[n]->getId());
			setNodeState(copyNode, nodeStates[0]);
		}
		else
		{
			setNodeState(copyNode, 0);
		}
	}
	return solution;
}


/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char** argv)
{
  cout << "************************************************************************************************" << endl;
  cout << "*       TraitRELAX program for detecting phenotype related changed in selective pressure       *" << endl;
  cout << "************************************************************************************************" << endl;
  cout << endl;

  try
  { 
    /* process input from params file */
    BppApplication bppml(args, argv, "bppml");

    /* process character data */
    const BinaryAlphabet* balpha = new BinaryAlphabet();
    VectorSiteContainer* charData = processCharacterData(&bppml, balpha);

    /* process tree */
    Tree* tree = processTree(&bppml);
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(tree))->getNodes();

	/* compute the likelihood function */
	TransitionModel* characterModel = setCharacterModel(&bppml, charData, balpha, tree);
	DiscreteDistribution *rDist = new ConstantRateDistribution();
	DRHomogeneousTreeLikelihood * characterTreeLikelihood = new DRHomogeneousTreeLikelihood (*tree, *charData, characterModel, rDist);

    /* compute ASR ML history */
    Tree* solution = computeASRMLSolution(characterTreeLikelihood);
    /* write solution */
    updateStatesInNodesNames(solution); // compute the likelihood given the mapping
    // write newick string to file
    string treeStr = TreeTools::treeToParenthesis(*solution);
    string filepath = ApplicationTools::getStringParameter("output.tree.file", bppml.getParams(), "", "", true, 1);
	ofstream file (filepath);
	file << treeStr << "\n";
	file.close();
	
    /* write the partition data to labels file */
	setPartition(&bppml, solution);

    delete balpha;
    delete charData;
	delete characterModel;
	delete rDist;
	delete characterTreeLikelihood;
    delete tree;
	delete solution;

    bppml.done();

  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}