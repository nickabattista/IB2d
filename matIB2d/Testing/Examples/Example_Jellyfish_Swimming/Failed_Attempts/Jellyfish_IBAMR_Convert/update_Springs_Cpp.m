#include "update_springs_freq.h"
#include <ibamr/IBSpringForceSpec.h>
 
void
update_springs_freq(
           tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
           LDataManager* const l_data_manager,
           const double current_time,
           const double dt)
{
  static const double pi = 4*atan(1);
  const int finest_ln = hierarchy->getFinestLevelNumber();
  double freq = 0.5;
  double theta;
 
  theta=freq*current_time*2*pi;
 
  
  // Find out the Lagrangian index ranges. 
  const std::pair<int,int>& muscle_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);
     
  // Get the LMesh (which we assume to be associated with the finest level of
  // the patch hierarchy).  Note that we currently need to update both "local"
  // and "ghost" node data.
  Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
  vector<LNode*> nodes;
  nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
  nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());
 
  // Update the spring lengths in their associated spring specs.
  tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
  for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
      LNode* node_idx = *it;
      IBSpringForceSpec* spring_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
         
      if (spring_spec == NULL) continue;  // skip to next node
 
      // Here we update the resting length of the spring
      //
      // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
      //        resting_length      is the resting length of the current Lagrangian point that is the "master index'
      //                                Since there may be more than one spring associated with each point, it's a vector.
      //        resting_length[0]   is the resting length of the first spring
      //        resting_length[1]   would be the resting length of the second spring associated with that Lagrangian point.
      //
      // In this example, the resting length is increased by 0.01*dt each time step.
 
      const int lag_idx = node_idx->getLagrangianIndex();
      double& phase = spring_spec->getParameters()[0][1];
       
      //getParameters()[a][b];  a = which spring (1st on node, 2nd, etc.), b=0: stiffness, b=1: resting lengths
        
      //Note that you can also getStiffnesses
 
      if (muscle_lag_idxs.first <= lag_idx && lag_idx < muscle_lag_idxs.second)
    {
        
      phase=sin(theta);
      //resting_length = Resting_aph[s1] ;
 
      /*if (t0>.5){ 
        spring_stiffness=0;
      }else{
        spring_stiffness=5000;
             
        }*/
    }
      //          if((current_time-floor(current_time))==t0)
      //{
      //  resting_length[0]+=0;
      //}
    }
 
return;
}// update_springs