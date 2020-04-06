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
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOFrequenciesSetFormat.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/Likelihood/JointLikelihoodFunction.h>

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
    SiteContainerTools::changeGapsToUnknownCharacters(*sites); // convert gaps to unknown characters (as done in bppML.cpp)
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

TransitionModel* setCharacterModel(BppApplication* bppml, VectorSiteContainer* charData, const BinaryAlphabet* alphabet, Tree* tree, DRTreeParsimonyScore* mpData)
{
    // create the model
    SubstitutionModel* model = new TwoParameterBinarySubstitutionModel(alphabet);

    // compute the maximum parsimony score and set the lower and upper bounds on mu (=rate) as mp/tree_size, 2*mp/tree_size 
    VDouble treeBranches = dynamic_cast<TreeTemplate<Node>*>(tree)->getBranchLengths();
    double treeSize = 0;
    for (size_t i=0; i<treeBranches.size(); ++i)
    {
        treeSize += treeBranches[i];
    }
    double characterMuLb = mpData->getScore()/treeSize;
    double characterMuUb = 4*characterMuLb; // used 4*lb instead of 2*lb because factor of 2 is not enough (optimization converges to upper bound)
    
    // set the initial values of the model
    if (!ApplicationTools::getBooleanParameter("character_model.set_initial_parameters", bppml->getParams(), true, "", true, false))
    {
        // set the value of mu to be the middle of the interval
        model->setParameterValue(string("mu"),(characterMuLb + characterMuUb) / 2);
        dynamic_cast<TwoParameterBinarySubstitutionModel*>(model)->setMuBounds(characterMuLb, characterMuUb);
        // estimate the initial frequencies as observedPseudoCount with pseudocount as 1 to avoid possible case of frequency = 0
        model->setFreqFromData(dynamic_cast<const SequenceContainer&>(*charData), 1); // the second arguemnt stands for pesudocount 1
    }
    else
    {
      double mu = ApplicationTools::getDoubleParameter("character_model.mu", bppml->getParams(), 10);
      model->setParameterValue(string("mu"), mu);      
      double pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", bppml->getParams(), 0.5);
      map<int,double>frequencies;
      frequencies[0] = pi0;
      frequencies[1] = 1-pi0;
      model->setFreq(frequencies);
    }

    return dynamic_cast<TransitionModel*>(model);
    
}

/******************************************************************************/

void setMpPartition(BppApplication* bppml, DRTreeParsimonyScore* mpData, const VectorSiteContainer* characterData, TransitionModel* characterModel, Tree* tree)
{
  mpData->computeSolution();
  const Tree& solution = mpData->getTree();
  vector <const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>&>(solution)).getNodes();
  
  // set assignment to branches
  string character0NodesIds, character1NodesIds = "";
  for (size_t i=0; i<nodes.size(); ++i)
  {
    if (!tree->isRoot(static_cast<int>(i)))
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
  } 
  bppml->getParam("model1.nodes_id") = character0NodesIds.substr(0,character0NodesIds.length()-1);
  bppml->getParam("model2.nodes_id") = character1NodesIds.substr(0,character1NodesIds.length()-1);  
}

/******************************************************************************/

MixedSubstitutionModelSet* setSequenceModel(BppApplication* bppml, const VectorSiteContainer* codon_data, const CodonAlphabet* codonAlphabet, DRTreeParsimonyScore* mpData, const VectorSiteContainer* characterData, TransitionModel* characterModel, Tree* tree)
{
    if (!ApplicationTools::getBooleanParameter("sequence_model.set_initial_parameters", bppml->getParams(), true, "", true, false))
    {
      bppml->getParam("model1") = "RELAX(kappa=1,p=0.1,omega1=1.0,omega2=2.0,theta1=0.5,theta2=0.8,k=1,frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
      bppml->getParam("model2") = "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,k=1,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
    }
    bppml->getParam("nonhomogeneous")="general";
    bppml->getParam("nonhomogeneous.number_of_models") = "2";
    bppml->getParam("nonhomogeneous.stationarity") = "yes"; // constrain root frequencies to be the same as stationary (since RELAX is a time reversible model, this should not cause issues)
    
    // set likelihood computation to restrict the same selective regime for each site along the phylogeny
    bppml->getParam("site.number_of_paths") = "2";                               // the 3rd path mapping omega3 in the branches under chatacter states 0 and 1 is imlies the the other two paths
    bppml->getParam("site.path1") = "model1[YN98.omega_1]&model2[YN98.omega_1]"; // map omega1 in the branches under character state 0 (=model1) to omega1 in the branches under character state 1 (=model2) 
    bppml->getParam("site.path2") = "model1[YN98.omega_2]&model2[YN98.omega_2]"; // do the same for omega2
    
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml->getParams(), "Standard", "", true, true);
    GeneticCode* gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);
    
    // set initial partition, based on maximum parsimony
    setMpPartition(bppml, mpData, characterData, characterModel, tree); // the partition is set on tree
    
    // create the set of models
    MixedSubstitutionModelSet* modelSet = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(codonAlphabet, gCode, codon_data, bppml->getParams()));
    return modelSet;
}



/******************************************************************************/


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
    BppApplication traitRELAX(args, argv, "traitRELAX");

    /* process character data */
    const BinaryAlphabet* balpha = new BinaryAlphabet();
    VectorSiteContainer* charData = processCharacterData(&traitRELAX, balpha);

    /* process tree */
    Tree* tree = processTree(&traitRELAX);
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(tree))->getNodes();

    /* process codon alignment */
    const CodonAlphabet* calpha = getCodonAlphabet();
    VectorSiteContainer* seqData = processCodonAlignment(&traitRELAX, calpha);

    /* compute the maximum parsimony  for the purpose of setting bounds on the rate parameter of the character model and an intial tree partition for the starting point */
    DRTreeParsimonyScore* mpData  = new DRTreeParsimonyScore(*tree,  dynamic_cast<const SiteContainer&>(*charData)); 

    /* set the character model */
    TransitionModel* charModel = setCharacterModel(&traitRELAX, charData, balpha, tree, mpData);

    /* set the joint likelihood function instance */
    DiscreteDistribution* rDist = new ConstantRateDistribution();

    /* set the joint likelihood function */
    ParameterList emptyParametersList;

    /* for stochastic mappings */
    cout << "**** Analysis based on 1000 mappings ****" << endl;
    cout << "** Sampling based analysis **" << endl;
    traitRELAX.getParam("character.num_of_mappings")  ="1000";
	traitRELAX.getParam("character.use_analytic_mapping")  ="0";
    MixedSubstitutionModelSet* seqModel = setSequenceModel(&traitRELAX, seqData, calpha, mpData, charData, charModel, tree);
    JointLikelihoodFunction* traitRELAXLikelihoodFunction_0 = new JointLikelihoodFunction(&traitRELAX, tree, charData, charModel, seqData, seqModel, rDist, true);
    traitRELAXLikelihoodFunction_0->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(0)); // only compute likelihood given the gimulated parameters
    traitRELAXLikelihoodFunction_0->setHypothesis(JointLikelihoodFunction::Hypothesis(1));
    traitRELAXLikelihoodFunction_0->fireParameterChanged(emptyParametersList);

    cout << "***************************************" << endl;
    cout << " ** Analytic rewards based analysis ** " << endl;
    traitRELAX.getParam("character.use_analytic_mapping")  ="1";
    traitRELAXLikelihoodFunction_0->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(0)); // only compute likelihood given the gimulated parameters
    traitRELAXLikelihoodFunction_0->setHypothesis(JointLikelihoodFunction::Hypothesis(1));
    traitRELAXLikelihoodFunction_0->fireParameterChanged(emptyParametersList);
    traitRELAX.getParam("character.use_analytic_mapping")  ="0";
    cout << "***************************************" << endl;

    // compute likelihood given the true history and write to file in debug dir
    string trueHistoryPath = ApplicationTools::getStringParameter("true_history.tree.file", traitRELAX.getParams(), "", "", true, 1);
    traitRELAX.getParam("input.tree.file") = trueHistoryPath;
    Tree* trueHistory = processTree(&traitRELAX);
    string bgLabels = ApplicationTools::getStringParameter("true_history.model1.nodes_id", traitRELAX.getParams(), "", "", true, 1);
    traitRELAX.getParam("model1.nodes_id") = bgLabels;
    string fgLabels = ApplicationTools::getStringParameter("true_history.model2.nodes_id", traitRELAX.getParams(), "", "", true, 1);
    traitRELAX.getParam("model2.nodes_id") = fgLabels;
    TransitionModel* newCharModel = setCharacterModel(&traitRELAX, charData, balpha, trueHistory, mpData);
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", traitRELAX.getParams(), "Standard", "", true, true);
    GeneticCode* gCode = SequenceApplicationTools::getGeneticCode(calpha->getNucleicAlphabet(), codeDesc);
    MixedSubstitutionModelSet* newSeqModel = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(calpha, gCode, seqData, traitRELAX.getParams()));  
    JointLikelihoodFunction* RELAXLikelihoodFunction = new JointLikelihoodFunction(&traitRELAX, trueHistory, charData, newCharModel, seqData, newSeqModel, rDist);
    RELAXLikelihoodFunction->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(0)); // only compute likelihood given the simulated parameters
    RELAXLikelihoodFunction->setHypothesis(JointLikelihoodFunction::Hypothesis(1));
    const RNonHomogeneousMixedTreeLikelihood* seqTl = RELAXLikelihoodFunction->getSequenceLikelihoodFunction();
    const HomogeneousTreeLikelihood* charTl = RELAXLikelihoodFunction->getCharacterLikelihoodFunction();
    string debugDir = ApplicationTools::getStringParameter("output.debug.dir", traitRELAX.getParams(), "", "", true, 1); 
    string filePath = debugDir + "true_history_likelihood.txt";
    ofstream myfile (filePath);
    if (myfile.is_open())
    {
        myfile << -1*charTl->getValue() + -1*seqTl->getValue() << "\n";
    }
    myfile.close();
    ApplicationTools::displayResult("Joint model log likelihood: ", TextTools::toString(((-1*charTl->getValue()) + (-1*seqTl->getValue())), 15));
            
    // now, conduct stablity Analysis for the sampling based likelihood computation
    string filepPath = debugDir + "sampling_expected_history_approach_logl_stability.txt";
    ofstream myfie (filepPath);
    if (myfile.is_open())
    {
      for (size_t i=0; i<1000; ++i)
      {
        traitRELAXLikelihoodFunction_0->fireParameterChanged(emptyParametersList);
        myfile << -1*traitRELAXLikelihoodFunction_0->getValue() << "\n";
      }
    }

    delete seqModel;
    delete balpha;
    delete rDist;
    delete charData;
    delete charModel;
    delete newCharModel;
    delete mpData;

    delete calpha;
    delete seqData;
    delete seqModel;
    delete newSeqModel;
    
    delete traitRELAXLikelihoodFunction_0;
    //delete traitRELAXLikelihoodFunction_1;
    //delete traitRELAXLikelihoodFunction_2;
    delete RELAXLikelihoodFunction;
    delete tree;
    delete trueHistory;

    traitRELAX.done();

  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}