function dydt = funk(t,y)


dydt = [ y(2); (1/12)*(-4*y(2)-3*y(1)) ];

end