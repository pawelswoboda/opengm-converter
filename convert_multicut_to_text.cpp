#include <fstream>
#include "opengm/opengm.hxx"
#include "opengm/graphicalmodel/graphicalmodel.hxx"
#include "opengm/graphicalmodel/graphicalmodel_hdf5.hxx"
#include "opengm/operations/adder.hxx"
#include "assert.h"

int main(int argc, char** argv)
{
   assert(argc == 3); // first arg input in opengm format, second is output in text format
   const std::string input_file = argv[1];
   const std::string output_file = argv[2];

   typedef double                                                               ValueType;          // type used for values
   typedef size_t                                                               IndexType;          // type used for indexing nodes and factors (default : size_t)
   typedef size_t                                                               LabelType;          // type used for labels (default : size_t)
   typedef opengm::Adder                                                        OpType;             // operation used to combine terms
   typedef opengm::ExplicitFunction<ValueType,IndexType,LabelType>              ExplicitFunction;   // shortcut for explicite function
   typedef opengm::PottsFunction<ValueType,IndexType,LabelType>                 PottsFunction;      // Potts function
   typedef opengm::meta::TypeListGenerator<ExplicitFunction,PottsFunction>::type              FunctionTypeList;   // list of all function the model cal use (this trick avoids virtual methods) - here only one
   typedef opengm::DiscreteSpace<IndexType, LabelType>                          SpaceType;          // type used to define the feasible statespace
   typedef opengm::GraphicalModel<ValueType,OpType,FunctionTypeList,SpaceType>  Model;              // type of the model
   typedef Model::FunctionIdentifier                                            FunctionIdentifier; // type of the function identifier

   Model gm;
   opengm::hdf5::load(gm, input_file,"gm");

   std::ofstream file_stream(output_file, std::ofstream::out);
   file_stream << "MULTICUT\n";

   for(std::size_t e=0; e<gm.numberOfFactors(); ++e) {
     assert(gm[e].numberOfVariables() == 2);
     const std::size_t i = gm.variableOfFactor(e,0);
     const std::size_t j = gm.variableOfFactor(e,1);
     assert(i < j);

     const double cost = gm[e](std::array<std::size_t,2>({1,0}).begin()) - gm[e](std::array<std::size_t,2>({0,0}).begin());

     file_stream << i << " " << j << " " << cost << "\n";
   }

   file_stream.close();
}
