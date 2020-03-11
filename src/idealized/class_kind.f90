module class_kind
implicit none
integer,  parameter   :: pd = 15
integer,  parameter   :: rd = 307
integer,  parameter   :: dp = selected_real_kind(pd,rd), wp = dp
integer,  parameter   :: qp = selected_real_kind (32)
end module class_kind