MODULE APHchl_Module
   USE Numeric_Kinds_Module
   SAVE
   INTEGER(KIND=IW_Kind):: naph     ! The number of channels in the APH propagation region. 
   INTEGER(KIND=IW_Kind):: neigmin_sum  ! The number of channels in calculated APH propagation region.
   INTEGER(KIND=IW_Kind):: Neigmin_Save(0:1000)  ! The number of channels in calculated APH propagation region.
   NAMELIST/aphchl/ naph
END MODULE APHchl_Module
