// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>       /* pow */

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
    vector<Node*> nodes = ttree.getNodes();
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

const CodonAlphabet* getCodonAlphabet()
{
  std::map<std::string, std::string> params;
  params["alphabet"] = "Codon(letter=DNA)";
  params["genetic_code"] = "Standard";
  const Alphabet* alphabet = SequenceApplicationTools::getAlphabet(params, "", false);
  const CodonAlphabet* codonAlphabet = dynamic_cast<const CodonAlphabet*>(alphabet);
  return codonAlphabet;
}

// curently inialized with random parameters - in the future, initial parameters will be based on maximum parsimony and mayne initial optimizations in the case of the sequence model
SubstitutionModelSet* setSequenceModel(BppApplication* bppml, const VectorSiteContainer* codon_data, const CodonAlphabet* codonAlphabet, double omega, double k)
{
    map<string,string> params;
    params["model1"] = "YN98(kappa=1,omega=" + TextTools::toString(omega) + ",frequencies=F0)";
    // set the omega in model2 to be omega^k
    params["model2"] = "YN98(kappa=YN98.kappa_1,omega=" + TextTools::toString(pow(omega,k)) + ",frequencies=F0)";
    params["nonhomogeneous"] = "general";
    params["nonhomogeneous.number_of_models"] = "2";
    params["nonhomogeneous.stationarity"] = "yes"; // constrain root frequencies to be the same as stationary (since RELAX is a time reversible model, this should not cause issues)
    
    GeneticCode* gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), "Standard");
    
    // create the set of models
    SubstitutionModelSet* modelSet = dynamic_cast<SubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(codonAlphabet, gCode, codon_data, params));
    return modelSet;
}


void homogenize(vector<SubstitutionModel*>  models, vector<double>Probs, const CodonAlphabet* codonAlphabet)
{
    // step 0: compute synto and synfrom as - the minimal indices i,j of synonymous transition for which Q matrices of model 1 and model 2 have rate > 0 (that is, Q(i,j) > 0)
    // unforetunatly, the roginial data memebers are protected and thus must be comptated from scratch rather than extracted
    // emulated from YNGP_M2.cpp lines 110-121
    vector<int> supportedChars = models[0]->getAlphabetStates();
    GeneticCode* gc = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), "Standard");
    size_t synfrom = 0;
	size_t synto = 0;
    for (synfrom = 1; synfrom < supportedChars.size(); ++synfrom)
    {
        for (synto = 0; synto < synfrom; ++synto)
        {
        if (gc->areSynonymous(supportedChars[synfrom], supportedChars[synto])
            && (models[0]->Qij(synfrom, synto) != 0)
            && (models[1]->Qij(synfrom, synto) != 0))
            break;
        }
        if (synto < synfrom)
        break;
    }
    if (synto == supportedChars.size())
        throw Exception("Impossible to find synonymous codons");

    // step 1: compute the vector of new rates as 1 / Q(synfrom, synto) for each model
    // emulated from YNGP_M2.cpp lines 136-142
    Vdouble Rates;
    Rates.push_back(1 / models[0]->Qij(synfrom, synto));
    Rates.push_back(1 / models[1]->Qij(synfrom, synto));
    Rates.push_back(1 / models[2]->Qij(synfrom, synto)); 

    // step 2: normalize the rates from step 1 as if all the models belong to the same mixture model: compute the weighted sum of rates
    // emulated from AbstractMixedSubstitutionModel.cpp lines 231-244
    double sum = 0;
    double sP = 0;
    for (unsigned int i = 0; i < Rates.size(); i++)
    {
        sum += Rates[i] * Probs[i];
        sP += Probs[i];
    }
    sum /= sP;

    // set the normalized rate of each model as its new rate
    for (unsigned int i=0; i<Rates.size(); i++)
    {
        Rates[i] *= 1 / sum;
        models[i]->setRate(Rates[i]);
    }
}

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char** argv)
{
  cout << "************************************************************************************************" << endl;
  cout << "*                                     TraitRELAX simulator                                     *" << endl;
  cout << "************************************************************************************************" << endl;
  cout << endl;

  try
  { 

    /* process input from params file */
    BppApplication simulationParams(args, argv, "simulationParams"); // the parameters file should consists of the partition in the form "model1.nodes_id=...\nmodel2.nodes_id=..."

    // extract simulation parameters
    size_t numOfSites = ApplicationTools::getIntParameter("num_of_sites", simulationParams.getParams(), 300, "", true, 1);
    size_t numOfReplicates = ApplicationTools::getIntParameter("num_of_replicates", simulationParams.getParams(), 100, "", true, 1);

    // sequence model parameters (for now, using F0 for frequencies)
    // omegas and their weights corresponds to the model radical simulation in RELAX paper (see figure 4)
    vector<double>omegas;
    omegas.push_back(ApplicationTools::getDoubleParameter("omega0", simulationParams.getParams(), 0.1, "", true, 1));
    omegas.push_back(ApplicationTools::getDoubleParameter("omega1", simulationParams.getParams(), 1, "", true, 1));
    omegas.push_back(ApplicationTools::getDoubleParameter("omega2", simulationParams.getParams(), 2, "", true, 1));
    vector<double> omegaWeights;
    omegaWeights.push_back(ApplicationTools::getDoubleParameter("p0", simulationParams.getParams(), 0.5, "", true, 1));
    omegaWeights.push_back(ApplicationTools::getDoubleParameter("p1", simulationParams.getParams(), 0.4, "", true, 1));
    omegaWeights.push_back(1-omegaWeights[0]-omegaWeights[1]);
    // process k values to simulate with
    vector<double> k_values = ApplicationTools::getVectorParameter<double>("k_values", simulationParams.getParams(), ',', ':', "", "");
    

    // process the base tree
    Tree* ttree = processTree(&simulationParams);
    TreeTemplate<Node>* tree = dynamic_cast<TreeTemplate<Node>*>(ttree);
    giveNamesToInternalNodes(ttree);
    // write the base tree to the main analysis path 
    PhylogeneticsApplicationTools::writeTree(*tree, simulationParams.getParams());

    DiscreteDistribution* rdist = new ConstantRateDistribution();
    vector<string> seqNames = tree->getLeavesNames();

    // generate an empty dataset and use it to create RELAX
    const CodonAlphabet* codonAlphabet = getCodonAlphabet();
    VectorSiteContainer seqDataContainter(seqNames, codonAlphabet);

    // simulate a codon alignment according to RELAX model
    // RELAX is a mixture of 3 non homogenous models: one for each omega. Thus, to simulate under RELAX, we choose for each site the model to simulate it under, based on multimodel distrubution over the 3 models

    // create the 3 branch models
    double k;
    for (unsigned int i=0; i < k_values.size(); ++i)
    {  
        k = k_values[i];
        for (unsigned int rep=0; rep < numOfReplicates; ++rep)
        {
            SubstitutionModelSet* seqModel_1 = setSequenceModel(&simulationParams, &seqDataContainter, codonAlphabet, omegas[0], k);
            SubstitutionModelSet* seqModel_2 = setSequenceModel(&simulationParams, &seqDataContainter, codonAlphabet, omegas[1], k);
            SubstitutionModelSet* seqModel_3 = setSequenceModel(&simulationParams, &seqDataContainter, codonAlphabet, omegas[2], k);

            /* emulate homogenization */
            vector<SubstitutionModel*> bgModels;
            bgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_1->getModel(0)));
            bgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_2->getModel(0)));
            bgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_3->getModel(0)));
            homogenize(bgModels, omegaWeights, codonAlphabet);
            vector<SubstitutionModel*> fgModels;
            fgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_1->getModel(1)));
            fgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_2->getModel(1)));
            fgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_3->getModel(1)));
            homogenize(fgModels, omegaWeights, codonAlphabet);

            // simulate sites
            NonHomogeneousSequenceSimulator* seqSimulator_1 = new NonHomogeneousSequenceSimulator(seqModel_1, rdist, tree);
            NonHomogeneousSequenceSimulator* seqSimulator_2 = new NonHomogeneousSequenceSimulator(seqModel_2, rdist, tree);
            NonHomogeneousSequenceSimulator* seqSimulator_3 = new NonHomogeneousSequenceSimulator(seqModel_3, rdist, tree);
            VectorSiteContainer seqData(seqNames, codonAlphabet);
			vector<size_t> omegaOrders;   // holds 1,2,3

            Site* seqSite;
            
            for (size_t site=0; site<numOfSites; ++site)
            {          
            
			omegaOrders.clear();
            omegaOrders.push_back(1);
            omegaOrders.push_back(2);
            omegaOrders.push_back(3);
			
            omegaWeights.clear();
            omegaWeights.push_back(0.5);
            omegaWeights.push_back(0.4);
            omegaWeights.push_back(1-omegaWeights[0]-omegaWeights[1]);
			
            size_t modelOfSite = RandomTools::pickOne(omegaOrders, omegaWeights);
            if (modelOfSite == 1)
            {
                seqSite = seqSimulator_1->simulateSite();
            }
            else if (modelOfSite == 2)
            {
                seqSite = seqSimulator_2->simulateSite();
            }
            else
            {
                seqSite = seqSimulator_3->simulateSite();
            }
            seqSite->setPosition(static_cast<int>(site));
            seqData.addSite(*seqSite, false);
            }

            // delete the models
            delete seqModel_1;
            delete seqModel_2;
            delete seqModel_3;
			
			omegaWeights.clear();
            omegaWeights.push_back(0.5);
            omegaWeights.push_back(0.4);
            omegaWeights.push_back(1-omegaWeights[0]-omegaWeights[1]);

            // write the simulated sequence data
            SequenceApplicationTools::writeSequenceFile(*(dynamic_cast<SequenceContainer*>(&seqData)), simulationParams.getParams(), "", false, 1);
        }
    }
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}