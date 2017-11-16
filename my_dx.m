function out = my_dx(in)
  out = .50*(circshift(in,[ 1, 0]) - circshift(in,[-1, 0]));