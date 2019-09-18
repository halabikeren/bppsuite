// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
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

const CodonAlphabet *getCodonAlphabet()
{
  map<string, string> alphabetSettings;
  alphabetSettings["alphabet"] = "Codon(letter=DNA)";
  const Alphabet *alphabet = SequenceApplicationTools::getAlphabet(alphabetSettings, "", false);
  const CodonAlphabet *codonAlphabet = dynamic_cast<const CodonAlphabet *>(alphabet);
  return codonAlphabet;
}

/******************************************************************************/

VectorSiteContainer *processCharacterData(BppApplication *bppml, const BinaryAlphabet *alphabet)
{
  string charDataFilePath = ApplicationTools::getAFilePath("input.character.file", bppml->getParams(), true, true, "", true, "none", 1);
  string sequenceFormat = ApplicationTools::getStringParameter("input.character.format", bppml->getParams(), "Fasta()", "", true, 1);
  BppOAlignmentReaderFormat bppoReader(1);
  unique_ptr<IAlignment> iAln(bppoReader.read(sequenceFormat));
  map<string, string> args(bppoReader.getUnparsedArguments());
  ApplicationTools::displayResult("character data file ", charDataFilePath);
  ApplicationTools::displayResult("chatacter data format ", iAln->getFormatName());
  const SequenceContainer *charCont = iAln->readAlignment(charDataFilePath, alphabet);
  VectorSiteContainer *sites = new VectorSiteContainer(*dynamic_cast<const OrderedSequenceContainer *>(charCont));
  delete charCont;
  return sites;
}

/******************************************************************************/

VectorSiteContainer *processCodonAlignment(BppApplication *bppml, const CodonAlphabet *codonAlphabet)
{
  VectorSiteContainer *allSites = SequenceApplicationTools::getSiteContainer(codonAlphabet, bppml->getParams());            // here, gaps will be converted to unknown characters
  VectorSiteContainer *sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml->getParams(), "", true, false); // convert the alignemnt to condensed format of unique sites
  delete allSites;                                                                                                          // delete the non-condenced intance of the sequence data
  SiteContainerTools::changeGapsToUnknownCharacters(*sites);                                                                // convert gaps to unknown characters (as done in bppML.cpp)
  return sites;
}

/******************************************************************************/

void giveNamesToInternalNodes(Tree *tree)
{
  TreeTemplate<Node> *ttree = dynamic_cast<TreeTemplate<Node> *>(tree);
  vector<Node *> nodes = ttree->getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (!nodes[i]->hasName())
      nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
  }
}

/******************************************************************************/

Tree *processTree(BppApplication *bppml)
{
  Tree *tree = PhylogeneticsApplicationTools::getTree(bppml->getParams());
  giveNamesToInternalNodes(tree);

  string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml->getParams(), "Input", "", true, 1); // process tree branches lengths
  string cmdName;
  map<string, string> cmdArgs;
  KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs); // this line process cmdName - which dictates weather the tree should be processed as is, or ultrameterized
  if (cmdName == "Clock")                                         // if needed, turn the tree into an ultrametric tree
  {
    TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
  }
  return tree;
}

/******************************************************************************/

TransitionModel *setCharacterModel(BppApplication *bppml, VectorSiteContainer *charData, const BinaryAlphabet *alphabet, Tree *tree, DRTreeParsimonyScore *mpData)
{
  // create the model
  // TO DO: ADD HERE PROCESSING OF INITIAL CHARACTER MODEL PARAMETERS AND PADD TO CONSTRUCTOR
  // extract the user initial value of k for potential later use
  double init_mu = ApplicationTools::getDoubleParameter("character_model.mu", bppml->getParams(), 1);
  double init_pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", bppml->getParams(), 0.5);
  SubstitutionModel *model = new TwoParameterBinarySubstitutionModel(alphabet, init_mu, init_pi0);

  // compute the maximum parsimony score and set the lower and upper bounds on mu (=rate) as mp/tree_size, 2*mp/tree_size
  VDouble treeBranches = dynamic_cast<TreeTemplate<Node> *>(tree)->getBranchLengths();
  double treeSize = 0;
  for (size_t i = 0; i < treeBranches.size(); ++i)
  {
    treeSize += treeBranches[i];
  }
  double characterMuLb = mpData->getScore() / treeSize;
  double characterMuUb = 4 * characterMuLb; // used 4*lb instead of 2*lb because factor of 2 is not enough (optimization converges to upper bound)

  // set the initial values of the model
  if (!ApplicationTools::getBooleanParameter("character_model.set_initial_parameters", bppml->getParams(), true, "", true, false))
  {
    // set the value of mu to be the middle of the interval
    model->setParameterValue(string("mu"), (characterMuLb + characterMuUb) / 2);
    dynamic_cast<TwoParameterBinarySubstitutionModel *>(model)->setMuBounds(characterMuLb, characterMuUb);
    // estimate the initial frequencies as observedPseudoCount with pseudocount as 1 to avoid possible case of frequency = 0
    model->setFreqFromData(dynamic_cast<const SequenceContainer &>(*charData), 1); // the second arguemnt stands for pesudocount 1
  }
  else
  {
    double mu = ApplicationTools::getDoubleParameter("character_model.mu", bppml->getParams(), 10);
    model->setParameterValue(string("mu"), mu);
    if (mu < characterMuLb)
    {
      dynamic_cast<TwoParameterBinarySubstitutionModel *>(model)->setMuBounds(mu - 0.001, characterMuUb);
    }
    else if (mu > characterMuUb)
    {
      dynamic_cast<TwoParameterBinarySubstitutionModel *>(model)->setMuBounds(characterMuLb, mu + 0.001);
    }
    double pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", bppml->getParams(), 0.5);
    map<int, double> frequencies;
    frequencies[0] = pi0;
    frequencies[1] = 1 - pi0;
    model->setFreq(frequencies);
  }

  return dynamic_cast<TransitionModel *>(model);
}

/******************************************************************************/

void setMpPartition(BppApplication *bppml, DRTreeParsimonyScore *mpData, const VectorSiteContainer *characterData, TransitionModel *characterModel, Tree *tree)
{
  mpData->computeSolution();
  const Tree &solution = mpData->getTree();
  vector<const Node *> nodes = (dynamic_cast<const TreeTemplate<Node> &>(solution)).getNodes();

  // set assignment to branches
  string character0NodesIds, character1NodesIds = "";
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (!tree->isRoot(static_cast<int>(i)))
    {
      int nodeState = dynamic_cast<const BppInteger *>(nodes[i]->getNodeProperty("state"))->getValue();
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
  bppml->getParam("model1.nodes_id") = character0NodesIds.substr(0, character0NodesIds.length() - 1);
  bppml->getParam("model2.nodes_id") = character1NodesIds.substr(0, character1NodesIds.length() - 1);
}

/******************************************************************************/

MixedSubstitutionModelSet *setSequenceModel(BppApplication *bppml, const VectorSiteContainer *codon_data, const CodonAlphabet *codonAlphabet, DRTreeParsimonyScore *mpData, const VectorSiteContainer *characterData, TransitionModel *characterModel, Tree *tree)
{
  if (!ApplicationTools::getBooleanParameter("sequence_model.set_initial_parameters", bppml->getParams(), true, "", true, false))
  {
    bppml->getParam("model1") = "RELAX(kappa=1,p=0.1,omega1=1.0,omega2=2.0,theta1=0.5,theta2=0.8,k=1,frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
    bppml->getParam("model2") = "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,k=1,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
  }
  bppml->getParam("nonhomogeneous") = "general";
  bppml->getParam("nonhomogeneous.number_of_models") = "2";
  bppml->getParam("nonhomogeneous.stationarity") = "yes"; // constrain root frequencies to be the same as stationary (since RELAX is a time reversible model, this should not cause issues)

  // set likelihood computation to restrict the same selective regime for each site along the phylogeny
  bppml->getParam("site.number_of_paths") = "2";                               // the 3rd path mapping omega3 in the branches under chatacter states 0 and 1 is imlies the the other two paths
  bppml->getParam("site.path1") = "model1[YN98.omega_1]&model2[YN98.omega_1]"; // map omega1 in the branches under character state 0 (=model1) to omega1 in the branches under character state 1 (=model2)
  bppml->getParam("site.path2") = "model1[YN98.omega_2]&model2[YN98.omega_2]"; // do the same for omega2

  string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml->getParams(), "Standard", "", true, true);
  GeneticCode *gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);

  // set initial partition, based on maximum parsimony
  setMpPartition(bppml, mpData, characterData, characterModel, tree); // the partition is set on tree

  // create the set of models
  MixedSubstitutionModelSet *modelSet = dynamic_cast<MixedSubstitutionModelSet *>(PhylogeneticsApplicationTools::getSubstitutionModelSet(codonAlphabet, gCode, codon_data, bppml->getParams()));
  return modelSet;
}

/******************************************************************************/

void initCharacterGrid(VVDouble &grid, double pi0Lb, double pi0Ub, double muLb, double muUb, uint gridSize)
{
  double stepSize;
  double next;

  /* define a set of possible assignments for rate */
  if (muLb < 0 && muUb > 1) // in case of a range that allows both relative character substitution rate lower and greater than 1, set the jumps in the grid to be uneven in order to sample same amount of rates < 1 and rates > 1
  {
    double logUb = log(muUb);
    double logLb = log(muLb);
    stepSize = (logUb - logLb) / (gridSize - 1);
    next = logLb;
    for (size_t i = 0; i < gridSize; ++i)
    {
      grid[0][i] = exp(next);
      next += stepSize;
    }
  }
  else
  {
    stepSize = (muUb - muLb) / (gridSize - 1);
    next = muLb;
    for (size_t j = 0; j < gridSize; ++j)
    {
      grid[0][j] = next;
      next += stepSize;
    }
  }

  /* define a set of possible assignments for kappa */
  stepSize = (pi0Ub - pi0Lb) / (gridSize - 1);
  next = pi0Lb;
  for (size_t k = 0; k < gridSize; ++k)
  {
    grid[1][k] = next;
    next += stepSize;
  }
}

/******************************************************************************/

VVDouble getNeighbors(VVDouble &grid, double currentMu, double currentPi0)
{
  VVDouble neighbors;
  neighbors.clear();
  neighbors.resize(2, VDouble(3));
  // get the location of the current point in the grid
  uint muPos = 0;
  uint pi0Pos = 0;
  for (uint i = 0; i < grid.size(); i++)
  {
    if (grid[0][i] == currentMu)
    {
      muPos = i;
      break;
    }
  }
  for (uint j = 0; j < grid.size(); j++)
  {
    if (grid[1][j] == currentPi0)
    {
      pi0Pos = j;
      break;
    }
  }
  // return a double-list of all its surronding points
  int muMin, muMax, pi0Min, pi0Max;
  muMin = muPos - 1;
  muMax = muPos + 1;
  if (muPos == 0)
  {
    muMin = muMin;
  }
  if (muPos == grid.size())
  {
    muMax = muPos;
  }
  pi0Min = pi0Pos - 1;
  pi0Max = pi0Pos + 1;
  if (pi0Pos == 0)
  {
    pi0Min = pi0Min;
  }
  if (pi0Pos == grid.size())
  {
    pi0Max = pi0Pos;
  }
  for (uint i = muMin; i < static_cast<uint>(muMax); i++)
  {
    for (uint j = pi0Min; j < static_cast<uint>(pi0Max); j++)
    {
      neighbors[0][i - muMin] = grid[0][i];
      neighbors[1][j - pi0Min] = grid[1][j];
    }
  }
  return neighbors;
}

/******************************************************************************/

void optimizeAlternativeCharacterModelByGrid(map<string, double> &optimalValues, JointLikelihoodFunction *jointlikelihoodFunction, TransitionModel *characterModel, uint verbose = 1, uint gridSize = 10, bool firstCycle = true)
{
  cout << "* Optimizng joint likelihood function with respect to character parameters using grid *\n"
       << endl;

  // set the grid bounds on the rate and kappa parameters of the character model
  const IntervalConstraint *muBounds = dynamic_cast<const IntervalConstraint *>(characterModel->getParameter("mu").getConstraint());
  double muLb = muBounds->getLowerBound() + 0.0001;
  double muUb = muBounds->getUpperBound() - 0.0001;
  double pi0Lb = 0.1;
  double pi0Ub = 0.9; // altered from 999 to avoid radical values (-> failure to simulate history along a branch) set the upper bound on kappa such that pi0 = 1/(kappa+1) in [0.001,0.999] -> kappa in [0,999] (exclude 0 and 1 from the bouds to avoid absorving states in the character model)
  VVDouble grid;
  grid.clear();
  grid.resize(2, VDouble(gridSize));
  initCharacterGrid(grid, pi0Lb, pi0Ub, muLb, muUb, gridSize);

  // scan the grid and for each assignment of (rate, kappa), compute the joint likelihood and update the pair yieldling the best likelihood accordingly
  Parameter mu = jointlikelihoodFunction->getParameter("TwoParameterBinary.mu");
  Parameter pi0 = jointlikelihoodFunction->getParameter("TwoParameterBinary.pi0");
  double currentMu = mu.getValue();
  double currentPi0 = pi0.getValue();
  double currentLogL = -1 * jointlikelihoodFunction->getValue();
  optimalValues["mu"] = currentMu;
  optimalValues["pi0"] = currentPi0;
  optimalValues["log_likelihood"] = currentLogL;
  if (firstCycle)
  {
    for (size_t i = 0; i < gridSize; ++i)
    {
      currentMu = grid[0][i];
      mu.setValue(currentMu); // set the value of the rate
      for (size_t j = 0; j < gridSize; ++j)
      {
        currentPi0 = grid[1][j];
        pi0.setValue(currentPi0); // set the value of kappa
        ParameterList paramsToUpdate;
        paramsToUpdate.addParameter(mu);
        paramsToUpdate.addParameter(pi0);
        jointlikelihoodFunction->setParametersValues(paramsToUpdate); // wrong values of parameters are set here. trigger likelihood computation
        currentLogL = -1 * jointlikelihoodFunction->getValue();
        if (currentLogL > optimalValues["log_likelihood"])
        {
          optimalValues["log_likelihood"] = currentLogL;
          optimalValues["mu"] = currentMu;
          optimalValues["pi0"] = currentPi0;
        }
      }
    }
  }
  else // i case of a non initial cycle, only traverse the neighbors of the current point in the likelihood surface
  {
    VVDouble neighbors = getNeighbors(grid, currentMu, currentPi0);
    for (uint i = 0; i < neighbors.size(); i++)
    {
      currentMu = grid[0][i];
      mu.setValue(currentMu); // set the value of the rate
      currentPi0 = grid[1][i];
      pi0.setValue(currentPi0); // set the value of kappa
      ParameterList paramsToUpdate;
      paramsToUpdate.addParameter(mu);
      paramsToUpdate.addParameter(pi0);
      jointlikelihoodFunction->setParametersValues(paramsToUpdate); // wrong values of parameters are set here. trigger likelihood computation
      currentLogL = -1 * jointlikelihoodFunction->getValue();
      if (currentLogL > optimalValues["log_likelihood"])
      {
        optimalValues["log_likelihood"] = currentLogL;
        optimalValues["mu"] = currentMu;
        optimalValues["pi0"] = currentPi0;
      }
    }
  }
}

/******************************************************************************/

// 27.8.18, note to future: use powell for two dimentional brent
void optimizeAlternativeCharacterModelByBrent(map<string, double> &optimalValues, JointLikelihoodFunction *jointlikelihoodFunction, TransitionModel *characterModel, uint verbose = 1)
{
  cout << "* Optimizing joint likelihood function with respect to character parameters using one dimentional brent *\n"
       << endl;

  // // set the brent one dimontional optimizer
  // optimalValues["log_likelihood"] = -INFINITY;
  BrentOneDimension *characterParametersOptimizer = new BrentOneDimension(jointlikelihoodFunction);
  characterParametersOptimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);
  characterParametersOptimizer->getStopCondition()->setTolerance(0.01);               // set the tolerance to be slighly less strict to account for the instability of the joint likelihood function
  characterParametersOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO); // 3.4.19 - keren - OTHERWISE, CONSTRAINT IS EXCEEDED
  characterParametersOptimizer->setProfiler(0);
  characterParametersOptimizer->setMessageHandler(0);
  characterParametersOptimizer->setVerbose(1);

  // optimize the joint model with respect to pi0
  ParameterList pi0;
  pi0.addParameter(jointlikelihoodFunction->getParameter("TwoParameterBinary.pi0"));
  const IntervalConstraint *pi0Bounds = dynamic_cast<const IntervalConstraint *>(jointlikelihoodFunction->getParameter("TwoParameterBinary.pi0").getConstraint());
  characterParametersOptimizer->setInitialInterval(pi0Bounds->getLowerBound(), pi0Bounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
  characterParametersOptimizer->init(pi0);
  characterParametersOptimizer->optimize();

  optimalValues["pi0"] = jointlikelihoodFunction->getParameter("TwoParameterBinary.pi0").getValue();

  // optimize the model with respect to the mu
  ParameterList mu;
  mu.addParameter(jointlikelihoodFunction->getParameter("TwoParameterBinary.mu"));
  const IntervalConstraint *muBounds = dynamic_cast<const IntervalConstraint *>(jointlikelihoodFunction->getParameter("TwoParameterBinary.mu").getConstraint());
  characterParametersOptimizer->setInitialInterval(muBounds->getLowerBound() + 0.0001, muBounds->getUpperBound() - 0.0001);
  characterParametersOptimizer->init(mu);
  characterParametersOptimizer->optimize();
  optimalValues["mu"] = jointlikelihoodFunction->getParameter("TwoParameterBinary.mu").getValue();
  optimalValues["log_likelihood"] = jointlikelihoodFunction->getCharacterLikelihoodFunction()->getValue();
  ApplicationTools::displayResult("Character log Likelihood after brent: ", TextTools::toString(-1 * jointlikelihoodFunction->getCharacterLikelihoodFunction()->getValue(), 15));
  delete characterParametersOptimizer;
}

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char **argv)
{
  cout << "************************************************************************************************" << endl;
  cout << "*       TraitRELAX program for detecting phenotype related changed in selective pressure       *" << endl;
  cout << "************************************************************************************************" << endl;
  cout << endl;

  try
  {

    /* process input from params file */
    BppApplication traitRELAX(args, argv, "traitRELAX");
    uint verbose = static_cast<unsigned int>(ApplicationTools::getIntParameter("verbose", traitRELAX.getParams(), 1));

    /* process character data */
    const BinaryAlphabet *balpha = new BinaryAlphabet();
    VectorSiteContainer *charData = processCharacterData(&traitRELAX, balpha);

    /* process tree */
    Tree *tree = processTree(&traitRELAX);
    vector<Node *> nodes = (dynamic_cast<TreeTemplate<Node> *>(tree))->getNodes();

    /* process codon alignment */
    const CodonAlphabet *calpha = getCodonAlphabet();
    VectorSiteContainer *seqData = processCodonAlignment(&traitRELAX, calpha);

    /* compute the maximum parsimony  for the purpose of setting bounds on the rate parameter of the character model and an intial tree partition for the starting point */
    DRTreeParsimonyScore *mpData = new DRTreeParsimonyScore(*tree, dynamic_cast<const SiteContainer &>(*charData));

    /* set the character model */
    TransitionModel *charModel = setCharacterModel(&traitRELAX, charData, balpha, tree, mpData);

    /* set the sequence model */
    MixedSubstitutionModelSet *seqModel = setSequenceModel(&traitRELAX, seqData, calpha, mpData, charData, charModel, tree);

    /* set the joint likelihood function instance */
    DiscreteDistribution *rDist = new ConstantRateDistribution();

    JointLikelihoodFunction *traitRELAXLikelihoodFunction = new JointLikelihoodFunction(&traitRELAX, tree, charData, charModel, seqData, seqModel, rDist);
    cout << "\n**** Initial parameters ****" << endl;
    map<string, double> userInitialValues = traitRELAXLikelihoodFunction->getModelParameters();
    ParameterList emptyParametersList;

    /* fit the null model: separate optimization of the character model and the sequence model, where the selection intensity parameter k is 1 */
    cout << "\n**** Null model fitting ****" << endl;
    traitRELAX.startTimer();
    traitRELAXLikelihoodFunction->setHypothesis(JointLikelihoodFunction::Hypothesis(0));
    traitRELAXLikelihoodFunction->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(3));
    traitRELAXLikelihoodFunction->fireParameterChanged(emptyParametersList);
    traitRELAX.done();

    /* report the log likelihood of the null model */
    cout << "\n**** Null model after optimization ****\n"
         << endl;
    traitRELAXLikelihoodFunction->getModelParameters();
    vector<double> nullLoglBySite = traitRELAXLikelihoodFunction->getLikelihoodForEachSite();

    /* fit the alternative model: sequencial optimization of the character model and then sequence model, given an expected history based on the character model */
    cout << "\n**** Alternative model fitting ****" << endl;
    traitRELAXLikelihoodFunction->setHypothesis(JointLikelihoodFunction::Hypothesis(1));
    traitRELAXLikelihoodFunction->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(2));

    // set starting point by optimizing sequence paTotal execution time:rameters why relaxing the constraint on the selection intensity parameter, so that the maximum parsmony based partiton will be considered
    cout << "\n** Setting starting point by relaxing constraint on the selection intensity parameter for model2 and using maximum parsimony based partition **\n"
         << endl;
    traitRELAX.startTimer();
    // optimize with starting point 1: null inference result
    cout << "Starting point 1: null fitting result" << endl;
    RNonHomogeneousMixedTreeLikelihood *sequenceTl = traitRELAXLikelihoodFunction->getSequenceLikelihoodFunction();
    OptimizationTools::optimizeTreeScale(sequenceTl, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
    PhylogeneticsApplicationTools::optimizeParameters(sequenceTl, sequenceTl->getParameters(), traitRELAX.getParams());
    map<string, double> nullSpResult = traitRELAXLikelihoodFunction->getModelParameters();
    // optimize with starting point 2: user input values
    cout << "\nStarting point 2: user input values" << endl;
    traitRELAXLikelihoodFunction->scaleSequenceTree(1);
    sequenceTl = traitRELAXLikelihoodFunction->getSequenceLikelihoodFunction();
    for (map<string, double>::iterator it = userInitialValues.begin(); it != userInitialValues.end(); it++)
    {
      if ((it->first.find("RELAX") != std::string::npos))
      {
        traitRELAXLikelihoodFunction->setParameterValue(it->first, it->second);
      }
    }
    OptimizationTools::optimizeTreeScale(sequenceTl, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
    PhylogeneticsApplicationTools::optimizeParameters(sequenceTl, sequenceTl->getParameters(), traitRELAX.getParams());
    map<string, double> userSpResult = traitRELAXLikelihoodFunction->getModelParameters(true);
    // set the initial values according to the winning starting point
    string winningSp = "user initial values";
    map<string, double> winningSpValues = userSpResult;
    if (-1.0 * userSpResult["Overall Log likelihood"] < -1.0 * nullSpResult["Overall Log likelihood"])
    {
      winningSp = "null fitting result";
      winningSpValues = nullSpResult;
    }
    cout << "\n\nWinning starting point: " << winningSp << endl;
    for (map<string, double>::iterator it = winningSpValues.begin(); it != winningSpValues.end(); it++)
    {
      if ((it->first.find("RELAX") != std::string::npos))
      {
        traitRELAXLikelihoodFunction->setParameterValue(it->first, it->second);
      }
      else if (it->first.compare("sequenceScalingFactor") == 0)
      {
        traitRELAXLikelihoodFunction->scaleSequenceTree(it->second);
      }
    }
    // compute the initial log likelihood
    cout << "* Starting iterative two-layer optimzation of the alternative model *" << endl;
    traitRELAXLikelihoodFunction->getSequenceLikelihoodFunction()->computeTreeLikelihood();
    // report final values and complete step
    traitRELAXLikelihoodFunction->getModelParameters(true);
    traitRELAX.done();
    traitRELAXLikelihoodFunction->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(0)); // no optimization of sequece parameters is required at these stage
    double currentAlternativeOverallLogL = -1.0 * traitRELAXLikelihoodFunction->getValue();
    double previousAlternativeOverallLogL = currentAlternativeOverallLogL;
    map<string, double> bestModelParameters = traitRELAXLikelihoodFunction->getModelParameters(false);
    double bestLogl = currentAlternativeOverallLogL;
    double optimizationCyclesNum = 0;
    double tolerance = 0.1;
    do
    {
      traitRELAX.startTimer();
      /* optimize the character model while fixing the sequence model with respect to the joint model */
      cout << "\n** Step 1: fix sequence model parameters, optimize character model parameters **\n"
           << endl;
      traitRELAXLikelihoodFunction->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(0)); // no optimization of sequece parameters is required at these stage
      string characterOptimizationMethod = ApplicationTools::getStringParameter("optimization.character_method", traitRELAX.getParams(), "brent");
      map<string, double> optimalCharacterParameters;
      optimalCharacterParameters.clear();
      if (characterOptimizationMethod.compare("grid") == 0)
      {
        uint gridSize = static_cast<unsigned int>(ApplicationTools::getIntParameter("optimization.grid_size", traitRELAX.getParams(), 10));
        optimizeAlternativeCharacterModelByGrid(optimalCharacterParameters, traitRELAXLikelihoodFunction, charModel, verbose, gridSize);
      }
      else // use brent
      {
        optimizeAlternativeCharacterModelByBrent(optimalCharacterParameters, traitRELAXLikelihoodFunction, charModel, verbose);
      }

      /* report optimal character parameters */
      cout << "\n** Optimal character parameters: **\n"
           << endl;
      ApplicationTools::displayResult("Mu", TextTools::toString(optimalCharacterParameters["mu"]));
      ApplicationTools::displayResult("Pi0", TextTools::toString(optimalCharacterParameters["pi0"]));
      traitRELAX.done();

      /* optimize the sequence model while fixing the character model with respect to the joint model */
      cout << "\n** Step 2: fix character model parameters, optimize sequence model parameters **\n"
           << endl;
      traitRELAX.startTimer();
      traitRELAXLikelihoodFunction->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(2)); // trigger optimization of the model
      traitRELAXLikelihoodFunction->fireParameterChanged(emptyParametersList);
      map<string, double> currentmodelParameters = traitRELAXLikelihoodFunction->getModelParameters();
      previousAlternativeOverallLogL = currentAlternativeOverallLogL;
      currentAlternativeOverallLogL = -1 * traitRELAXLikelihoodFunction->getValue();
      traitRELAX.done();
      if (currentAlternativeOverallLogL > bestLogl)
      {
        bestModelParameters = currentmodelParameters;
        bestLogl = -1.0 * bestModelParameters["Overall Log likelihood"];
      }
      ApplicationTools::displayResult("Optimization cycle", TextTools::toString(optimizationCyclesNum + 1));
      ApplicationTools::displayResult("Difference between current and previous joint log likelihood", TextTools::toString(currentAlternativeOverallLogL - previousAlternativeOverallLogL));
      optimizationCyclesNum = optimizationCyclesNum + 1;
      cout << "\n\n************************************************************************************************\n\n"
           << endl;
    } while (currentAlternativeOverallLogL - previousAlternativeOverallLogL > tolerance);

    /* report the optimal log likelihood and parameters of the alternative model */
    cout << "\n**** Alternative model likelihood after optimization ****\n"
         << endl;
    for (map<string, double>::iterator it = bestModelParameters.begin(); it != bestModelParameters.end(); it++)
    {
      if ((it->first.find("Log likelihood") != std::string::npos))
      {
        ApplicationTools::displayResult(it->first, TextTools::toString(-1.0 * it->second));
      }
      else
      {
        ApplicationTools::displayResult(it->first, TextTools::toString(it->second));
      }
    }
    // for monitoring efficacy of TraitRELAX
    if (-1.0 * bestModelParameters["Overall Log likelihood"] == -1.0 * winningSpValues["Overall Log likelihood"])
    {
      cout << "Starting point (given MP history) yields the best likelihood" << endl;
    }
    ApplicationTools::displayResult("Number of optimization cycles", TextTools::toString(optimizationCyclesNum));

    vector<double> altLoglBySite = traitRELAXLikelihoodFunction->getLikelihoodForEachSite();

    cout << "Computing log likelihood ratio by site" << endl;
    double logLR;
    cout << "site\tlog(LR)" << endl;
    for (size_t s=0; s<altLoglBySite.size(); ++s)
    {
        logLR = altLoglBySite[s] / nullLoglBySite[s];
        cout <<  s << "\t" << logLR << endl;
    }

    // free parameters
    delete balpha;
    delete rDist;
    delete charData;
    delete charModel;
    delete mpData;

    delete calpha;
    delete seqData;
    delete seqModel;

    delete traitRELAXLikelihoodFunction;
    delete tree;

    traitRELAX.done();
  }
  catch (exception &e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}