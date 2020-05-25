function [fire_name,save_name,prefix] = fire_choice()

fire_choice = input_num('which fire? Patch: [0], Camp: [1], Cougar: [3]',1);
cycle = input_num('Which cycle? ',0)
if fire_choice == 1
    fire_name = 'Camp fire';
    save_name = 'camp';
    prefix='../campTIFs/';
elseif fire_choice == 0
    fire_name = 'Patch Springs fire';
    save_name = 'patch';
    prefix='../TIFs/';
else
    fire_name = 'Cougar Creek fire';
    save_name = 'cougar';
    prefix = '../cougarTIFs/';
end

end