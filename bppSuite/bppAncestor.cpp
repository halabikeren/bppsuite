//
// File: bppAncestor.cpp
// Created by: Julien Dutheil
// Created on: Sep Wed 10 14:14 2008
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SiteTools.h>
#include <Seq/SequenceApplicationTools.h>
#include <Seq/SequenceContainerTools.h>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/DRHomogeneousTreeLikelihood.h>
#include <Phyl/DRNonHomogeneousTreeLikelihood.h>
#include <Phyl/PatternTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/MarginalAncestralStateReconstruction.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/RASTools.h>
#include <Phyl/Newick.h>
#include <Phyl/MarkovModulatedSubstitutionModel.h>
#include <Phyl/SubstitutionModelSet.h>
#include <Phyl/SubstitutionModelSetTools.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/AutoParameter.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  SequenceApplicationTools::printInputAlignmentHelp();
  PhylogeneticsApplicationTools::printInputTreeHelp();
  PhylogeneticsApplicationTools::printSubstitutionModelHelp();
  PhylogeneticsApplicationTools::printRateDistributionHelp();
  PhylogeneticsApplicationTools::printCovarionModelHelp();
  SequenceApplicationTools::printOutputSequenceHelp();
  *ApplicationTools::message << "asr.method                    | (marginal only for now)" << endl;
  *ApplicationTools::message << "asr.probabilites              | [bool] shall I print the posterior" << endl;
  *ApplicationTools::message << "                              | probabilities of each state?" << endl;
  *ApplicationTools::message << "asr.sample                    | [bool] shall I sample the posterior distribution" << endl;
  *ApplicationTools::message << "                              | instead of taking the max probability?" << endl;
  *ApplicationTools::message << "asr.sample.number             | [int>0] number of samples" << endl;
  *ApplicationTools::message << "output.file                   | file where to write site infos" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
}

void printParameters(const SubstitutionModel* model, ofstream& out)
{
  out << "model.name = " << model->getName() << endl;
  ParameterList pl = model->getParameters();
  for(unsigned int i = 0; i < pl.size(); i++)
    out << "model." << pl[i]->getName() << " = " << pl[i]->getValue() << endl;
}

void printParameters(const SubstitutionModelSet* modelSet, ofstream& out, map<string, string>& params)
{
  out << "nonhomogeneous = general" << endl;
  out << "nonhomogeneous.number_of_models = " << modelSet->getNumberOfModels() << endl;
  for(unsigned int i = 0; i < modelSet->getNumberOfModels(); i++)
  {
    const SubstitutionModel* model = modelSet->getModel(i);
    out << endl;
    out << "model" << (i+1) << ".name = " << model->getName() << endl;
    //ParameterList pl = model->getParameters();
    //for(unsigned int j = 0; j < pl.size(); j++)
    //  out << "model" << (i+1) << "." << pl[j]->getName() << " = " << pl[j]->getValue() << endl;
    vector<int> ids = modelSet->getNodesWithModel(i);
    out << "model" << (i+1) << ".nodes_id = " << ids[0];
    for(unsigned int j = 1; j < ids.size(); j++)
      out << "," << ids[j];
    out << endl;
    out << "model" << (i+1) << ".use_observed_freq = no" << endl;
  }
  ParameterList pl = modelSet->getParameters();
  ParameterList plroot = modelSet->getRootFrequenciesParameters();
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    if(plroot.getParameter(pl[i]->getName()) == NULL)
    {
      out << endl;
      string name = modelSet->getParameterModelName(pl[i]->getName());
      vector<unsigned int> models = modelSet->getModelsWithParameter(pl[i]->getName());
      out << "model" << (models[0] + 1) << "." << name << " = " << pl[i]->getValue() << endl;
      for(unsigned int j = 1; j < models.size(); j++)
      {
        out << "model" << (models[j] + 1) << "." << name << " = model" << (models[0] + 1) << "." << name << endl;
      }
    }
  }
 
  //Root frequencies:
  out << endl;
  out << "# Root frequencies:" << endl;
  if(plroot.size() == 1 && plroot[0]->getName() == "RootFreqtheta")
  {
    out << "nonhomogeneous.root_freq = initGC" << endl;
    out << "model.ancTheta = " << plroot[0]->getValue() << endl;
  }
  else
  {
    out << "nonhomogeneous.root_freq = init" << endl;
    vector<double> rootFreqs = modelSet->getRootFrequencies();
    for(unsigned int i = 0; i < rootFreqs.size(); i++)
      out << "model.anc" << modelSet->getAlphabet()->intToChar((int)i) << " = " << rootFreqs[i] << endl;
  }
}

void printParameters(const DiscreteDistribution* rDist, ofstream& out, map<string, string>& params)
{
  out << "rate_distribution = " << params["rate_distribution"] << endl;
  ParameterList pl = rDist->getParameters();
  for(unsigned int i = 0; i < pl.size(); i++)
    out << "rate_distribution." << pl[i]->getName() << " = " << pl[i]->getValue() << endl;
  bool invariant = (params["rate_distribution"].find("+invariant") != string::npos);
  out << "rate_distribution.classes_number = " << (rDist->getNumberOfCategories() + (invariant ? -1 : 0)) << endl;
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Ancestral Sequence Reconstruction, version 0.1.0     *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 10/09/08 *" << endl;
  cout << "*         B. Boussau                                             *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  ApplicationTools::startTimer();

  cout << "Parsing options:" << endl;
  
  map<string, string> params = AttributesTools::parseOptions(args, argv);

  Alphabet * alphabet = SequenceApplicationTools::getAlphabet(params, "", false);

  VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet, params);
  
  VectorSiteContainer * sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, params);
  delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
  
  // Get the initial tree
  TreeTemplate<Node> * tree = PhylogeneticsApplicationTools::getTree(params);
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  
  string treeWIdPath = ApplicationTools::getAFilePath("output.tree.path", params, false, false);
  if(treeWIdPath != "none")
  {
    vector<Node *> nodes = tree->getNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      if(nodes[i]->isLeaf())
        nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
      else
        nodes[i]->setBranchProperty("NodeId", String(TextTools::toString(nodes[i]->getId())));
    }
    Newick treeWriter;
    treeWriter.enableExtendedBootstrapProperty("NodeId");
    ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
    treeWriter.write(*tree, treeWIdPath);
    delete tree;
    cout << "BppML's done." << endl;
    exit(0);
  }

  DRTreeLikelihood *tl;
  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", params, "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  SubstitutionModel    * model    = NULL;
  SubstitutionModelSet * modelSet = NULL;
  DiscreteDistribution * rDist    = NULL;
  unsigned int nbStates;

  if(nhOpt == "no")
  {  
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, params);
    if(model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    }
    tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true);
    nbStates = model->getNumberOfStates();
  }
  else if(nhOpt == "one_per_branch")
  {
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, params);
    if(model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    }
    vector<double> rateFreqs;
    if(model->getNumberOfStates() != alphabet->getSize())
    {
      //Markov-Modulated Markov Model...
      unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                   // we should assume a rate distribution for the root also!!!  
    }
    FrequenciesSet * rootFreqs = PhylogeneticsApplicationTools::getFrequenciesSet(alphabet, sites, params, rateFreqs);
    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", params, ',', "");
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
    model = NULL;
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
    nbStates = model->getNumberOfStates();
  }
  else if(nhOpt == "general")
  {
    modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, sites, params);
    if(modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    }
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
    nbStates = modelSet->getNumberOfStates();
  }
  else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
  tl->initialize();
 
  delete tree;
    
  double logL = tl->getValue();
  if(isinf(logL))
  {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = tl->getBranchLengthsParameters();
    for(unsigned int i = 0; i < pl.size(); i++)
    {
      if(pl[i]->getValue() < 0.000001) pl[i]->setValue(0.000001);
    }
    tl->matchParametersValues(pl);
    logL = tl->getValue();
  }
  if(isinf(logL))
  {
    ApplicationTools::displayError("!!! Unexpected likelihood == 0.");
    ApplicationTools::displayError("!!! Looking at each site:");
    for(unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
      *ApplicationTools::error << "Site " << sites->getSite(i)->getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i) << endl;
    }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
    exit(-1);
  }
  tree = new TreeTemplate<Node>(*tl->getTree());

  // Write parameters to screen:
  ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
  ParameterList parameters = tl->getSubstitutionModelParameters();
  for(unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i]->getName(), TextTools::toString(parameters[i]->getValue()));
  }
  parameters = tl->getRateDistributionParameters();
  for(unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i]->getName(), TextTools::toString(parameters[i]->getValue()));
  }

  // Getting posterior rate class distribution:
  DiscreteDistribution * prDist = RASTools::getPosteriorRateDistribution(* tl);
  ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
  if(ApplicationTools::message) prDist->print(*ApplicationTools::message);
  ApplicationTools::displayMessage("\n");
  delete prDist;

  // Reconstruct ancestral sequences:
  string reconstruction = ApplicationTools::getStringParameter("asr.method", params, "marginal", "", true, false);
  ApplicationTools::displayResult("Ancestral state reconstruction method", reconstruction);
  bool probs = false;

  AncestralStateReconstruction *asr = NULL;
  bool probMethod = false;
  if(reconstruction == "marginal")
  {
    asr = new MarginalAncestralStateReconstruction(*tl);
    probMethod = true;
  }
  else
    throw Exception("Unknown ancestral state reconstruction method: " + reconstruction);

  if(probMethod)
  {
    probs = ApplicationTools::getBooleanParameter("asr.probabilities", params, false, "", true, false);
    ApplicationTools::displayResult("Output probabilities", probs ? "yes" : "no");
  }

  // Write infos to file:
  string outputFile = ApplicationTools::getAFilePath("output.file", params, false, false);
  if(outputFile != "none")
  {
    ApplicationTools::displayResult("Output file", outputFile);
    ofstream out(outputFile.c_str(), ios::out);
    vector<Node *> nodes = tree->getInnerNodes();
    unsigned int nbNodes = nodes.size();
    
    // Get the rate class with maximum posterior probability:
    vector<unsigned int> classes = tl->getRateClassWithMaxPostProbOfEachSite();
    // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
    Vdouble rates = tl->getPosteriorRateOfEachSite();
    // Get the ancestral sequences:
    vector<Sequence *> sequences(nbNodes);
    vector<VVdouble *> probabilities(nbNodes);

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");
    colNames.push_back("rc");
    colNames.push_back("pr");
    for(unsigned int i = 0; i < nbNodes; i++)
    {
      Node *node = nodes[i];
      colNames.push_back("max." + TextTools::toString(node->getId()));
      if(probs)
      {
        probabilities[i] = new VVdouble();
        //The cast will have to be updated when more probabilistic method will be available:
        sequences[i] = dynamic_cast<MarginalAncestralStateReconstruction *>(asr)->getAncestralSequenceForNode(node->getId(), probabilities[i], false);

        for(unsigned int j = 0; j < nbStates; j++)
        {
          colNames.push_back("prob." + TextTools::toString(node->getId()) + "." + alphabet->intToChar((int)j));
        }
      }
      else
      {
        sequences[i] = asr->getAncestralSequenceForNode(node->getId());
      }
    }

    //Now fill the table:
    vector<string> row(colNames.size());
    DataTable * infos = new DataTable(colNames);
    
    for(unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
      double lnL = tl->getLogLikelihoodForASite(i);
      const Site * currentSite = sites->getSite(i);
      int currentSitePosition = currentSite->getPosition();
      int isCompl = (SiteTools::isComplete(* currentSite) ? 1 : 0);
      int isConst = (SiteTools::isConstant(* currentSite) ? 1 : 0);
      row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
      row[1] = TextTools::toString(isCompl);
      row[2] = TextTools::toString(isConst);
      row[3] = TextTools::toString(lnL);
      row[4] = TextTools::toString(classes[i]);
      row[5] = TextTools::toString(rates[i]);

      unsigned int k = 6;
      for(unsigned int j = 0; j < nbNodes; j++)
      {
        row[k] = sequences[j]->getChar(i);
        k++;
        if(probs)
        {
          for(unsigned int l = 0; l < nbStates; l++)
          {
            row[k] = TextTools::toString((*probabilities[j])[i][l]);
            k++;
          }
        }
      }

      infos->addRow(row);
    }

    DataTable::write(*infos, out, "\t");

    delete infos;
  }
 
  SequenceContainer *asSites = NULL;
  if(probMethod)
  {
    bool sample = ApplicationTools::getBooleanParameter("asr.sample", params, false, "", true, false);
    ApplicationTools::displayResult("Sample from posterior distribution", sample ? "yes" : "no");
    if(sample)
    {
      unsigned int nbSamples = ApplicationTools::getParameter<unsigned int>("asr.sample.number", params, 1, "", true, false);
      asSites = new AlignedSequenceContainer(alphabet);
      for(unsigned int i = 0; i < nbSamples; i++)
      {
        ApplicationTools::displayGauge(i, nbSamples-1, '=');
        SequenceContainer *sampleSites = dynamic_cast<MarginalAncestralStateReconstruction *>(asr)->getAncestralSequences(true);
        vector<string> names = sampleSites->getSequencesNames();
        for(unsigned int j = 0; j < names.size(); j++)
          names[j] += "_" + TextTools::toString(i+1);
        sampleSites->setSequencesNames(names, false);
        SequenceContainerTools::append<AlignedSequenceContainer>(* dynamic_cast<AlignedSequenceContainer *>(asSites), *sampleSites);
        delete sampleSites;
      }
      *ApplicationTools::message << endl;
    }
    else
    {
      asSites = asr->getAncestralSequences();
    }
  }
  else
  {
    asSites = asr->getAncestralSequences();
  }
  SequenceApplicationTools::writeSequenceFile(*asSites, params);
  delete asSites;

  delete asr;
  delete alphabet;
  delete sites;
  if(model)    delete model;
  if(modelSet) delete modelSet;
  delete rDist;
  delete tl;
  delete tree;
  cout << "BppAncestor's done. Bye." << endl;
  ApplicationTools::displayTime("Total execution time:");
 
  }
  catch(exception & e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

