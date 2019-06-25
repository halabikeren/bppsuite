// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
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
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOFrequenciesSetFormat.h>


using namespace bpp;
using namespace std;

#include <map>
#include <vector>
//#include <boost/math/distributions/chi_squared.hpp>

/******************************************************************************/
/**************************** Auxiliary functions *****************************/
/******************************************************************************/


void printModelParameters(DiscreteRatesAcrossSitesTreeLikelihood* tl, MixedSubstitutionModelSet* modelSet)
{
  ParameterList parameters;
  for (size_t m = 0; m < modelSet->getNumberOfModels(); ++m) {
    ApplicationTools::displayMessage("\nmodel " + TextTools::toString(m+1) + "\n");
    TransitionModel* model = modelSet->getModel(m);
    parameters = model->getParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
  }
  parameters = tl->getRateDistributionParameters();
  for (size_t i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
}

VectorSiteContainer* process_alignment(Alphabet* alphabet, BppApplication bppml)
{
    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams()); // here, gaps will be converted to unknown characters
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), "", true, false); // convert the alignemnt to condensed format of unique sites
    delete allSites; // delete the non-condenced intance of the sequence data
    SiteContainerTools::changeGapsToUnknownCharacters(*sites); // convert gaps to unknown characters (as done in bppML.cpp)
    return(sites);

}

Tree* process_tree(BppApplication bppml)
{
    Tree* tree = PhylogeneticsApplicationTools::getTree(bppml.getParams());
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
    // try to add father to Internal4 (parent of Leaf3) - no such function - need to use detailsSimulation instance

    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml.getParams(), "Input", "", true, 1); // process tree branches lengths
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs); // this line process cmdName - which dictates weather the tree should be processed as is, or ultrameterized
    if (cmdName == "Clock") // if needed, turn the tree into an ultrametric tree
    {
      TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
    }
    return(tree);
}

MixedSubstitutionModelSet* setRELAXModel(BppApplication* bppml, const VectorSiteContainer* codon_data, const CodonAlphabet* codonAlphabet, Tree* tree)
{
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml->getParams(), "Standard", "", true, true);
    GeneticCode* gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);
    
    // create the set of models
    MixedSubstitutionModelSet* modelSet = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(codonAlphabet, gCode, codon_data, bppml->getParams()));
    return modelSet;
}


/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Test non-homogenous models    " << BPP_VERSION << "      *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  { 

    // process input from params file
    BppApplication bppml(args, argv, "bppML");

    map<string,string> parans = bppml.getParams(); // debug

    /* process alphabet type */
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }

    /* process alignment */
    VectorSiteContainer* sites = process_alignment(alphabet, bppml);

    /* process tree */
    Tree* tree = process_tree(bppml);
  
    /* set a branch site model */
    MixedSubstitutionModelSet* model = setRELAXModel(&bppml, sites, codonAlphabet, tree);

    /* generate a likelihood function */
    DiscreteDistribution* rDist = new ConstantRateDistribution();
    bppml.startTimer();
    RNonHomogeneousMixedTreeLikelihood* tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, true);
    tl->initialize();

    /*  save the inital value of k and then set RELAX.k_2 to 1 */
    // extract the user initial value of k for potential later use
    string FGModelInitialValues = ApplicationTools::getStringParameter("model2", bppml.getParams(), "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1)", "", true, true);
    string modelName = "RELAX";
    map<string, string> cmdArgs2;
    KeyvalTools::parseProcedure(FGModelInitialValues, modelName, cmdArgs2);
    double userInitialValue = TextTools::toDouble(cmdArgs2["k"]);
    tl->setParameterValue("RELAX.k_2", 1);

    /* compute likelihood */
    cout << "\nComputing intial log likelihood" << endl;
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));
    printModelParameters(tl, model); // debug - print model parameters
    bppml.done();

    /* fit the null */
    cout << "\nFitting the null model" << endl;
    bppml.startTimer();
    bppml.getParam("optimization.ignore_parameters") = "BrLen,RELAX.k_1,RELAX.k_2";
    PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams());
    double nullLogl = -tl->getValue();
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(nullLogl, 15));
    printModelParameters(tl, model); // debug - print model parameters
    bppml.done();

    /* fit the alternative */
    // set the initial value of k
    model->getModel(1)->setParameterValue("k", userInitialValue);
    cout << "\nFitting the alternative model" << endl;
    bppml.startTimer();
    bppml.getParam("optimization.ignore_parameters") = "BrLen,RELAX.k_1,RELAX.1_Full.theta_1,RELAX.1_Full.theta1_1,RELAX.1_Full.theta2_1,RELAX.2_Full.theta_1,RELAX.2_Full.theta1_1,RELAX.2_Full.theta2_1,RELAX.3_Full.theta_1,RELAX.3_Full.theta1_1,RELAX.3_Full.theta2_1"; // ignore frequency parameters to reduce optimization duration - results in one unit of ll reduction in optimality and 1 minutre reduction in duration
    PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams());
    double alternativeLogl = -tl->getValue();
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(alternativeLogl, 15));
    printModelParameters(tl, model); // debug - print model parameters
    bppml.done();

    /* do LRT */
    //doLRT(&bppml, nullLogl, alternativeLogl, 1);

    /* free resources */
    delete alphabet;
    delete sites;
    delete tree;
    if (model) delete model;
    delete rDist;
    delete tl;
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}