function out = my_dy(in)
  out = .50*(circshift(in,[ 0, 1]) - circshift(in,[0, -1]));