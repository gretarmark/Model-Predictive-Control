function Ts = get_Ts(Tref)

% Define sampling interval Ts (default Tref/10):
if Tref == 0,
  Ts = 1;
else
  Ts = Tref/10;
end

end