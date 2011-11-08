!> Definitions of all structures used during the construction of the Xpoint grid

module grid_xpoint_data
  implicit none
  integer, parameter :: n_flux_max = 1024
  integer, parameter :: n_tht_max  = 2048
  integer, parameter :: n_pieces   = 3

  type type_strategic_points                      !< type definition for strategic points
    real*8	      :: RLeftCorn_LowerInnerLeg,  ZLeftCorn_LowerInnerLeg  !< LeftCorn_LowerInnerLeg
    real*8	      :: RRightCorn_LowerInnerLeg, ZRightCorn_LowerInnerLeg !< RightCorn_LowerInnerLeg
    real*8	      :: RLeftCorn_LowerOuterLeg,  ZLeftCorn_LowerOuterLeg  !< LeftCorn_LowerOuterLeg  
    real*8	      :: RRightCorn_LowerOuterLeg, ZRightCorn_LowerOuterLeg !< RightCorn_LowerOuterLeg
    real*8	      :: RStrike_LowerInnerLeg,    ZStrike_LowerInnerLeg    !< Strike_LowerInnerLeg     
    real*8	      :: RStrike_LowerOuterLeg,    ZStrike_LowerOuterLeg    !< Strike_LowerOuterLeg    
    real*8	      :: RLeftCorn_UpperInnerLeg,  ZLeftCorn_UpperInnerLeg  !< LeftCorn_UpperInnerLeg 
    real*8	      :: RRightCorn_UpperInnerLeg, ZRightCorn_UpperInnerLeg !< RightCorn_UpperInnerLeg
    real*8	      :: RLeftCorn_UpperOuterLeg,  ZLeftCorn_UpperOuterLeg  !< LeftCorn_UpperOuterLeg  
    real*8	      :: RRightCorn_UpperOuterLeg, ZRightCorn_UpperOuterLeg !< RightCorn_UpperOuterLeg
    real*8	      :: RStrike_UpperInnerLeg,    ZStrike_UpperInnerLeg    !< Strike_UpperInnerLeg   
    real*8	      :: RStrike_UpperOuterLeg,    ZStrike_UpperOuterLeg    !< Strike_UpperOuterLeg   
    real*8	      :: RSecondStrike_InnerLeg,   ZSecondStrike_InnerLeg   !< SecondStrike_InnerLeg  
    real*8	      :: RSecondStrike_OuterLeg,   ZSecondStrike_OuterLeg   !< SecondStrike_OuterLeg  
    real*8	      :: RMiddle_LowerPrivate,     ZMiddle_LowerPrivate     !< Middle_LowerPrivate    
    real*8	      :: RMiddle_UpperPrivate,     ZMiddle_UpperPrivate     !< Middle_UpperPrivate     
    real*8	      :: RLimit_LowerOuterLeg,     ZLimit_LowerOuterLeg     !< Limit_LowerOuterLeg    
    real*8	      :: RLimit_LowerInnerLeg,     ZLimit_LowerInnerLeg     !< Limit_LowerInnerLeg    
    real*8	      :: RLimit_UpperOuterLeg,     ZLimit_UpperOuterLeg     !< Limit_UpperOuterLeg    
    real*8	      :: RLimit_UpperInnerLeg,     ZLimit_UpperInnerLeg     !< Limit_UpperInnerLeg 
    real*8	      :: angle_LowerLeft,          angle_LowerRight         !< Angle of horizontal line at lower Xpoint 
    real*8	      :: angle_UpperLeft,          angle_UpperRight         !< Angle of horizontal line at upper Xpoint
  end type type_strategic_points
  
  type type_new_points                      !< type definition for new grid points
    real*8 	      :: R_sep(n_tht_max), Z_sep(n_tht_max)
    real*8 	      :: R_max(n_tht_max), Z_max(n_tht_max)
    real*8 	      :: R_min(n_tht_max), Z_min(n_tht_max)
    real*8 	      :: RR_new(n_flux_max,n_tht_max),ZZ_new(n_flux_max,n_tht_max)
    real*8 	      :: s_flux(n_flux_max,n_tht_max),t_flux(n_flux_max,n_tht_max),t_tht(n_flux_max,n_tht_max)
    real*8 	      :: R_polar(n_pieces,4,n_tht_max),Z_polar(n_pieces,4,n_tht_max)
    integer           :: ielm_flux(n_flux_max,n_tht_max), k_cross(n_flux_max,n_tht_max)
  end type type_new_points
  
end module grid_xpoint_data

