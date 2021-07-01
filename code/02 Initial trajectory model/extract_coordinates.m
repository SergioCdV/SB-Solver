function  [rho, theta, z] = extract_coordinates(C, r0, tfapp)

rho = extract_coordinates_aux(C{1},r0,tfapp);
theta = extract_coordinates_aux(C{2},1,tfapp);
z = extract_coordinates_aux(C{3},r0,tfapp);

end

function coord = extract_coordinates_aux(C,adimdis, tfapp)

coord.o = C(:,1).*adimdis;
coord.D = C(:,2).*adimdis/tfapp;
coord.DD = C(:,3).*adimdis/tfapp^2;

end