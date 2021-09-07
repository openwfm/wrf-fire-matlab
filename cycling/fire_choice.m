function [fire_name,save_name,prefix,perim] = fire_choice(fire_choice)

if ~exist('fire_choice','var')
    fire_choice = input_num('which fire? Patch: [0], Camp: [1], Cougar: [3] other: [4]  Cal_2020 [5]',3,0);
end
%cycle = input_num('Which cycle? ',0,1)
if fire_choice == 1
    fire_name = 'Camp fire';
    save_name = 'camp';
    prefix='../campTIFs/';
    perim = '../PERIMs/camp_perims/';
elseif fire_choice == 0
    fire_name = 'Patch Springs fire';
    save_name = 'patch';
    prefix='../TIFs/';
    perim = '../PERIMs/patch_perims/';
elseif fire_choice == 4
    fire_name = 'CA fire';
    save_name = 'CA_2020';
    prefix='../creek/';
    perim = '../PERIMs/creek/'
elseif fire_choice == 5
    fire_name = 'CA fire';
    save_name = 'CA_2020';
    prefix='../calTIFS/';
    perim = '../PERIMs/august/'
else
    fire_name = 'Cougar Creek fire';
    save_name = 'cougar';
    prefix = '../cougarTIFs/';
    perim = '../PERIMs/cougar_perims/';
end

end