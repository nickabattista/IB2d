#include "update_target_point_positions.h"

// IBAMR INCLUDES
#include <ibamr/IBTargetPointForceSpec.h>

// IBTK INCLUDES
//#include <ibtk/LNodeIndexData.h>

// NAMESPACE
//using namespace IBAMR;

void
update_target_point_positions(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const lag_manager,
    const double current_time,
    const double dt)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // The velocity of the coral (dm/s).
    //static const double V = 0.1;
    
    //frames per second in data
    static const double F = 30.0; 
    //pixels per dm in data 
    static const double S = 1733.333333333333333333;
    //time in frames (mod one pulse)
    const double t =  F*fmod(current_time,2.30);
    const double tp = F*2.25; 
    const double tend = F*2.30; 
    //initial offset from center of the domain in dm
    static const double offset = 0.019; 
    static const double yoffset = -0.25; 

    static const double qc11 = -8.018203211249667e-13;
    static const double qc12 = 1.706310971452060e-10;
    static const double qc13 = -1.371526883228119e-8;
    static const double qc14 = 5.116506115621489e-7;
    static const double qc15 = -8.613944795521634e-6;
    static const double qc16 = 4.951095471423292e-5;
    static const double qc17 = -1.237254488643602e-4;
    static const double qc21 = -1.093988700839764e-10;
    static const double qc22 = 2.370737151797367e-8;
    static const double qc23 = -1.946684756263760e-6;
    static const double qc24 = 7.447999061932255e-5;
    static const double qc25 = -1.285757196003425e-3;
    static const double qc26 = 7.184871148722453e-3;
    static const double qc27 = -2.357284935085283e-2;    
    static const double qc31 = -3.111491615328405e-9;
    static const double qc32 = 6.922752279229525e-7; 
    static const double qc33 = -5.951659336196739e-5; 
    static const double qc34 = 2.461541908755476e-3; 
    static const double qc35 = -4.830999119980220e-2; 
    static const double qc36 = 3.293063238297988e-1;
    static const double qc37 = -1.893567487439901;
    static const double qc41 = 2.028150369106308e-8;
    static const double qc42 = -4.188949936239104e-6; 
    static const double qc43 = 3.237380275082298e-4;
    static const double qc44 = -1.158763275070524e-2; 
    static const double qc45 = 1.908244836752949e-1;
    static const double qc46 = -1.142372536532416; 
    static const double qc47 = 3.584246435337805;

    double qpoly1 = 0;
    double qpoly2 = 0;
    double qpoly3 = 0;
    double qpoly4 = 0; 

    if (t<=tp)
      {
	qpoly1=qc17+qc16*t+qc15*t*t+qc14*t*t*t+qc13*t*t*t*t+qc12*t*t*t*t*t+qc11*t*t*t*t*t*t; 
	qpoly2=qc27+qc26*t+qc25*t*t+qc24*t*t*t+qc23*t*t*t*t+qc22*t*t*t*t*t+qc21*t*t*t*t*t*t;
	qpoly3=qc37+qc36*t+qc35*t*t+qc34*t*t*t+qc33*t*t*t*t+qc32*t*t*t*t*t+qc31*t*t*t*t*t*t;
	qpoly4=qc47+qc46*t+qc45*t*t+qc44*t*t*t+qc43*t*t*t*t+qc42*t*t*t*t*t+qc41*t*t*t*t*t*t;
	//qpoly(1,:) = qcoeffs(1,7)+qcoeffs(1,6)*t+qcoeffs(1,5)*t.^2+qcoeffs(1,4)*t.^3+qcoeffs(1,3)*t.^4+qcoeffs(1,2)*t.^5+qcoeffs(1,1)*t.^6;
	//qpoly(2,:) = qcoeffs(2,7)+qcoeffs(2,6)*t+qcoeffs(2,5)*t.^2+qcoeffs(2,4)*t.^3+qcoeffs(2,3)*t.^4+qcoeffs(2,2)*t.^5+qcoeffs(2,1)*t.^6;
	//qpoly(3,:) = qcoeffs(3,7)+qcoeffs(3,6)*t+qcoeffs(3,5)*t.^2+qcoeffs(3,4)*t.^3+qcoeffs(3,3)*t.^4+qcoeffs(3,2)*t.^5+qcoeffs(3,1)*t.^6;
	//qpoly(4,:) = qcoeffs(4,7)+qcoeffs(4,6)*t+qcoeffs(4,5)*t.^2+qcoeffs(4,4)*t.^3+qcoeffs(4,3)*t.^4+qcoeffs(4,2)*t.^5+qcoeffs(4,1)*t.^6;
      }
    else
      {
	const double qpolytf1=qc17+qc16*tp+qc15*tp*tp+qc14*tp*tp*tp+qc13*tp*tp*tp*tp+qc12*tp*tp*tp*tp*tp+qc11*tp*tp*tp*tp*tp*tp;
	const double qpolytf2=qc27+qc26*tp+qc25*tp*tp+qc24*tp*tp*tp+qc23*tp*tp*tp*tp+qc22*tp*tp*tp*tp*tp+qc21*tp*tp*tp*tp*tp*tp;
	const double qpolytf3=qc37+qc36*tp+qc35*tp*tp+qc34*tp*tp*tp+qc33*tp*tp*tp*tp+qc32*tp*tp*tp*tp*tp+qc31*tp*tp*tp*tp*tp*tp;
	const double qpolytf4=qc47+qc46*tp+qc45*tp*tp+qc44*tp*tp*tp+qc43*tp*tp*tp*tp+qc42*tp*tp*tp*tp*tp+qc41*tp*tp*tp*tp*tp*tp;

	qpoly1 = qpolytf1*((tend-tp)-(t-tp))/(tend-tp)+qc17*(t-tp)/(tend-tp); 
	qpoly2 = qpolytf2*((tend-tp)-(t-tp))/(tend-tp)+qc27*(t-tp)/(tend-tp);
	qpoly3 = qpolytf3*((tend-tp)-(t-tp))/(tend-tp)+qc37*(t-tp)/(tend-tp);
	qpoly4 = qpolytf4*((tend-tp)-(t-tp))/(tend-tp)+qc47*(t-tp)/(tend-tp);
      }

    // Find out the Lagrangian index ranges.
    const std::pair<int,int>& coral2d_left_lag_idxs = lag_manager->getLagrangianStructureIndexRange(0, finest_ln);
    const std::pair<int,int>& coral2d_rght_lag_idxs = lag_manager->getLagrangianStructureIndexRange(1, finest_ln);

    // Get the patch data descriptor index for the LNodeIndexData.
    //const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Get the LMesh (which we assume to be associated with the finest level of                                                                  
    // the patch hierarchy).  Note that we currently need to update both "local"                                                                 
    // and "ghost" node data.                                                                                                                    
    Pointer<LMesh> mesh = lag_manager->getLMesh(finest_ln);
    vector<LNode*> nodes;
    nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
    nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

    // Update the target point positions in their associated target point force
    // specs.
    //tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
    //for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    //{
    //    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
    //    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
    //    const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
    //    for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
    //         it != idx_data->lnode_index_end(); ++it)
    //    {
    //        const LNodeIndex& node_idx = *it;
    //        const int lag_idx = node_idx.getLagrangianIndex();
    //        SAMRAI::tbox::Pointer<IBTargetPointForceSpec> force_spec = node_idx.getStashData<IBTargetPointForceSpec>();
    //        if (!force_spec.isNull())
    
    // Update the target point positions in their associated target point force                                                                  
    // specs.                                                                                                                                    
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
    for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
      {
        LNode* node_idx = *it;
        IBTargetPointForceSpec* force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
        if (force_spec == NULL) continue;  // skip to next node 
        //{
        // Here we update the position of the target point.
        //
        // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
        //        X_target     is the target position of the target point
        //        X_target[0]  is the x component of the target position
        //        X_target[1]  is the y component of the target position
        //
        // The target position is shifted to the left or right by the
        // increment dt*V
         
	//std::vector<double>& X_target = force_spec->getTargetPointPosition();

    const int lag_idx = node_idx->getLagrangianIndex();
	Point& X_target = force_spec->getTargetPointPosition();
	const double shifted_Xt1 = X_target[1]-yoffset; 

	X_target[0] = offset+qpoly4/S-qpoly3*shifted_Xt1+qpoly2*S*shifted_Xt1*shifted_Xt1-qpoly1*S*S*shifted_Xt1*shifted_Xt1*shifted_Xt1;

        if (coral2d_left_lag_idxs.first <= lag_idx && lag_idx < coral2d_left_lag_idxs.second)
        {	  
	  X_target[0] = -X_target[0]; 
	}

        if (coral2d_rght_lag_idxs.first <= lag_idx && lag_idx < coral2d_rght_lag_idxs.second)
        {
	  X_target[0] = X_target[0];
	}
      }     
    return;
}// update_target_point_positions

