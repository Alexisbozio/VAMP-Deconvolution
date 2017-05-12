function y=damp(x_new, x,damp_meas) 


y=damp_meas * x + (1 - damp_meas) * x_new;