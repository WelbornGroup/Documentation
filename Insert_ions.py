# Reference ion selection. Give reference ions' atom index here.
ref_1_index = 2
ref_2_index = 18
ref_3_index = 37
ref_4_index = 61

# Read the archive file.
with open ('ionchannel_comb.arc','r') as file:
    lines = file.readlines()

    # Collect the line numbver of each frame's strating point.
    start_of_frame = []
    for number, line in enumerate(lines):
        if "90.000000   90.000000   90.000000" in line:
            start_of_frame.append(number - 1)

    #Collect reference atoms coordinate from each frame.
    ref_1_selec = ref_1_index + 1
    ref_2_selec = ref_2_index + 1
    ref_3_selec = ref_3_index + 1
    ref_4_selec = ref_4_index + 1
    reference_1_coor = []
    reference_2_coor = []
    reference_3_coor = []
    reference_4_coor = []
    for i in range(len(start_of_frame)):
        current_ref_1 = []
        current_ref_2 = []
        current_ref_3 = []
        current_ref_4 = []
        for k in range(2,5):                                                                   # Get coordinates in current frame
            current_ref_1.append(float(lines[ref_1_selec + start_of_frame[i]].split()[k]))
            current_ref_2.append(float(lines[ref_2_selec + start_of_frame[i]].split()[k]))
            current_ref_3.append(float(lines[ref_3_selec + start_of_frame[i]].split()[k]))
            current_ref_4.append(float(lines[ref_4_selec + start_of_frame[i]].split()[k]))
        reference_1_coor.append(current_ref_1)
        reference_2_coor.append(current_ref_2)
        reference_3_coor.append(current_ref_3)
        reference_4_coor.append(current_ref_4)

    # Calculate mid point and collect into a list.
    x_mid_1 = []
    y_mid_1 = []
    z_mid_1 = []
    x_mid_2 = []
    y_mid_2 = []
    z_mid_2 = []   
    for i in range(len(reference_1_coor)):
        current_x_mid_1 = (reference_1_coor[i][0] + reference_2_coor[i][0]) / 2
        current_y_mid_1 = (reference_1_coor[i][1] + reference_2_coor[i][1]) / 2
        current_z_mid_1 = (reference_1_coor[i][2] + reference_2_coor[i][2]) / 2
        current_x_mid_2 = (reference_3_coor[i][0] + reference_4_coor[i][0]) / 2
        current_y_mid_2 = (reference_3_coor[i][1] + reference_4_coor[i][1]) / 2
        current_z_mid_2 = (reference_3_coor[i][2] + reference_4_coor[i][2]) / 2        
        # Set decimal point.
        current_x_round_1 = '%.6f'%current_x_mid_1
        current_y_round_1 = '%.6f'%current_y_mid_1
        current_z_round_1 = '%.6f'%current_z_mid_1
        current_x_round_2 = '%.6f'%current_x_mid_2
        current_y_round_2 = '%.6f'%current_y_mid_2
        current_z_round_2 = '%.6f'%current_z_mid_2        
        x_mid_1.append(current_x_round_1)
        y_mid_1.append(current_y_round_1)
        z_mid_1.append(current_z_round_1)
        x_mid_2.append(current_x_round_2)
        y_mid_2.append(current_y_round_2)
        z_mid_2.append(current_z_round_2)        

    # Atom info.
    atom_name = "Na+"
    atom_type = str(352)
    atom_index_1 = str(start_of_frame[1] - 1)
    atom_index_2 = str(start_of_frame[1])

    # Update header.
    tot_atoms = float(lines[0].split()[0])
    new_tot_atoms = str(int(tot_atoms + 2))
    header_rest = (" ",new_tot_atoms,"  ",lines[0].split()[1],"   ",lines[0].split()[2],"\n")
    header_comb = ''.join(header_rest)

    # Write to a new file.
    newfile = open ('ionchannel_NA.arc','w')
    for i in range(len(x_mid_1)):
        if float(x_mid_1[i]) < 0 and float(y_mid_1[i]) < 0 and float(z_mid_1[i]) < 0:
            insert_content_1 = (" ",atom_index_1,"  ",atom_name,"  ",x_mid_1[i],"  ",y_mid_1[i],"  ",z_mid_1[i],"   ",atom_type,"\n")
        elif float(x_mid_1[i]) > 0 and float(y_mid_1[i]) < 0 and float(z_mid_1[i]) < 0:
            insert_content_1 = (" ",atom_index_1,"  ",atom_name,"   ",x_mid_1[i],"  ",y_mid_1[i],"  ",z_mid_1[i],"   ",atom_type,"\n")
        elif float(x_mid_1[i]) < 0 and float(y_mid_1[i]) > 0 and float(z_mid_1[i]) < 0:
            insert_content_1 = (" ",atom_index_1,"  ",atom_name,"  ",x_mid_1[i],"   ",y_mid_1[i],"  ",z_mid_1[i],"   ",atom_type,"\n")
        elif float(x_mid_1[i]) < 0 and float(y_mid_1[i]) < 0 and float(z_mid_1[i]) > 0:
            insert_content_1 = (" ",atom_index_1,"  ",atom_name,"  ",x_mid_1[i],"  ",y_mid_1[i],"   ",z_mid_1[i],"   ",atom_type,"\n")
        elif float(x_mid_1[i]) > 0 and float(y_mid_1[i]) > 0 and float(z_mid_1[i]) < 0:
            insert_content_1 = (" ",atom_index_1,"  ",atom_name,"   ",x_mid_1[i],"   ",y_mid_1[i],"  ",z_mid_1[i],"   ",atom_type,"\n")
        elif float(x_mid_1[i]) > 0 and float(y_mid_1[i]) < 0 and float(z_mid_1[i]) > 0:
            insert_content_1 = (" ",atom_index_1,"  ",atom_name,"   ",x_mid_1[i],"  ",y_mid_1[i],"   ",z_mid_1[i],"   ",atom_type,"\n")
        elif float(x_mid_1[i]) < 0 and float(y_mid_1[i]) > 0 and float(z_mid_1[i]) > 0:
            insert_content_1 = (" ",atom_index_1,"  ",atom_name,"  ",x_mid_1[i],"   ",y_mid_1[i],"   ",z_mid_1[i],"   ",atom_type,"\n")
        else:
            insert_content_1 = (" ",atom_index_1,"  ",atom_name,"   ",x_mid_1[i],"   ",y_mid_1[i],"   ",z_mid_1[i],"   ",atom_type,"\n")

        if float(x_mid_2[i]) < 0 and float(y_mid_2[i]) < 0 and float(z_mid_2[i]) < 0:
            insert_content_2 = (" ",atom_index_2,"  ",atom_name,"  ",x_mid_2[i],"  ",y_mid_2[i],"  ",z_mid_2[i],"   ",atom_type,"\n")
        elif float(x_mid_2[i]) > 0 and float(y_mid_2[i]) < 0 and float(z_mid_2[i]) < 0:
            insert_content_2 = (" ",atom_index_2,"  ",atom_name,"   ",x_mid_2[i],"  ",y_mid_2[i],"  ",z_mid_2[i],"   ",atom_type,"\n")
        elif float(x_mid_2[i]) < 0 and float(y_mid_2[i]) > 0 and float(z_mid_2[i]) < 0:
            insert_content_2 = (" ",atom_index_2,"  ",atom_name,"  ",x_mid_2[i],"   ",y_mid_2[i],"  ",z_mid_2[i],"   ",atom_type,"\n")
        elif float(x_mid_2[i]) < 0 and float(y_mid_2[i]) < 0 and float(z_mid_2[i]) > 0:
            insert_content_2 = (" ",atom_index_2,"  ",atom_name,"  ",x_mid_2[i],"  ",y_mid_2[i],"   ",z_mid_2[i],"   ",atom_type,"\n")
        elif float(x_mid_2[i]) > 0 and float(y_mid_2[i]) > 0 and float(z_mid_2[i]) < 0:
            insert_content_2 = (" ",atom_index_2,"  ",atom_name,"   ",x_mid_2[i],"   ",y_mid_2[i],"  ",z_mid_2[i],"   ",atom_type,"\n")
        elif float(x_mid_2[i]) > 0 and float(y_mid_2[i]) < 0 and float(z_mid_2[i]) > 0:
            insert_content_2 = (" ",atom_index_2,"  ",atom_name,"   ",x_mid_2[i],"  ",y_mid_2[i],"   ",z_mid_2[i],"   ",atom_type,"\n")
        elif float(x_mid_2[i]) < 0 and float(y_mid_2[i]) > 0 and float(z_mid_2[i]) > 0:
            insert_content_2 = (" ",atom_index_2,"  ",atom_name,"  ",x_mid_2[i],"   ",y_mid_2[i],"   ",z_mid_2[i],"   ",atom_type,"\n")
        else:
            insert_content_2 = (" ",atom_index_2,"  ",atom_name,"   ",x_mid_2[i],"   ",y_mid_2[i],"   ",z_mid_2[i],"   ",atom_type,"\n")                
#        insert_content_1 = (" ",atom_index_1,"",atom_name," ",x_mid_1[i],"",y_mid_1[i],"",z_mid_1[i]," ",atom_type,"\n")
#        insert_content_2 = (" ",atom_index_2,"",atom_name," ",x_mid_2[i],"",y_mid_2[i],"",z_mid_2[i]," ",atom_type,"\n")
        ion_line_1 = "".join(insert_content_1)
        ion_line_2 = "".join(insert_content_2)
        frame_lines = "".join(lines[start_of_frame[i]+1:start_of_frame[i]+start_of_frame[1]])
        newfile.write(header_comb)
        newfile.write(frame_lines)
        newfile.write(ion_line_1)
        newfile.write(ion_line_2)
    newfile.close