MODULE FEM_Module
   USE Numeric_Kinds_Module
   SAVE
   ! All of the following parameters must be unique
   INTEGER, PARAMETER :: UniqueNode = -4    ! Unique Nodes       (center nodes)
   INTEGER, PARAMETER :: nonxst     = -3    ! Nonexistant nodes  
   INTEGER, PARAMETER :: cnstrain   = -2    ! Node constrained to another node
   INTEGER, PARAMETER :: free       =  0    ! Free nodes
   INTEGER, PARAMETER :: fixed      =  1    ! Fixed nodes (nodal value must be zero)
END MODULE FEM_Module