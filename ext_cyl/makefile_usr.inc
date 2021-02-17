
SRC_INC := ./src_inc
SRC_USR := ./src_usr
$(OBJDIR)/slp.o     :$(SRC_USR)/slp.f;                      $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
