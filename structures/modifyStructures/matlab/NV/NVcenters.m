%clear all;
%{
  This  function places nitrogen subtitutionals(ns), vacancies(v), vacancy pairs(v2), and nitrogen-vacancy(nv) Frenkel pairs
    in a carbon FCC diamond lattice. 
  This function requires one string and 6 INPUTS.
      
      name = the name of an an xyz file, minus the .xyz extension, in single quotes, e.g. '8x8x8'
      
      ns = the number of nitrogen defects randomly and uniformly disperesed in the lattice. The bond length between the nitrogen defect and each 1st
        nearest neighbor is stretched  to %28 longer than C-C bond according to 
        "J. Manuf. Mater. Process. 2017, 1, 6". 
        If only substitutionl defects are desired nv == 0.
        
      nv = the number of nv pairs in the lattice. A C nearest neighbor is deleted from a randomly, uniformly selected lattice site. 
        It is required that nv <= ns. The "bond" length between the vacancy
        and Carbon atoms is stretched to 1.63A and the dist. between the v
        and ns site is stretched to 1.69A according to "Phys. Rev. B 77,
        155206 (2008)."
      
      v = the number of vacany sites randomly and uniformly dispersed in the lattice. 
      
      v2 = the number of vacancy pairs to form from existing vacancy sites. It is required
      that v2 <= v.

      lammps = 0 or 1. 0 = DO NOT WRITE a lammps data file.
                       1 = write a lammps data file

      xyz = 0 or 1. 0 = DO NOT WRITE an xyz file.
                    1 = write an xyz file.
      
  ns centers are placed first, then nv sites, then v's, then v2's. By finding nn's to stretch N-C bonds, it is ensured
      that a v is not inadverently placed by a ns site.
%}

function NV(name, ns, nv, v, v2, lammps, xyz)

  d_nc = 0.25; %distance to shift C-N bond length from CC length
  d_nvc = 0.0486;
  d_nvn = 0.0833;
  CC_len = 1.5459; %CC bond length
  
  filename = strcat(name, '.xyz');
  
  fid = fopen(filename);
  
  sys = textscan(fid,'%d', 1); %ignore headers and find number of atoms
  sys = sys{1}; %number of atoms

  %%%%%%%%%%%%%%%%
  %Check inputs
  
  if nv > ns %nake sure nv <= ns
    except = 'nv must be less than or equal to ns';
    disp(except);
    return
  end
  
  if v2 > v %nake sure v <= v2
    except = 'v2 must be less than or equal to v';
    disp(except);
    return
  end
  
  flag = 0;

  if (ns < 0) 
    flag = 1;
  elseif (nv < 0)
    flag =1;
  elseif (v < 0)
    flag = 1;
  elseif (v2 < 0)
    flag = 1;
  elseif ((ns + nv + v + v2) > sys)
    flag = 1;
  end
  
  if flag == 1
    except = 'values must be non-negative and less than the total number of atoms';
    disp(except);
    return
  end
  
  %%%%%%%%%%%%%%%%%
  %Read data
  
  tmp = textscan(fid, '%f%d%d%d%f%d%d%d%f');  %ignore comments
  tmp = textscan(fid, '%2c%f%f%f', sys); %read in x, y, z coordinates

  fclose(fid);
  
  tmp_types = tmp{1,1}; %save atom types, not used
  len = length(tmp_types); %not used
  types = zeros(len,1); %not used 

  data = zeros(sys,5); %save ids, type, and x,y, and z coordinates to "data"
  ids = 1:1:sys;
  ids = transpose(ids);

  data(:,2) = 1;
  data(:,3) = tmp{1,2};
  data(:,4) = tmp{1,3};
  data(:,5) = tmp{1,4};

  [~,sort_ids] = sort(data(:,3)); %sort data by x coordinates
  data = data(sort_ids,:);
  data(:,1) = ids;
  
  %%%%%%%%%%%%%%%%%%%%%%
  %Place defects
  
  if ns ~= 0
    rand_ids = randperm(sys, ns); %find random ids and data for ns centers and remove those ids from sample population
    rand_ns = data(rand_ids,:); 
    data(rand_ids,:) = [];
    %%%pick a nearest neighbor. calculate the direction between them and tranlsate the ns atom to 28% longer than the C-C bond 
        
    if nv ~= 0 %remove nv sites from ns list
      rand_ids = randperm(size(rand_ns, 1), nv); %randomly choose some of the ns centers to form nv pairs
      rand_nv = rand_ns(rand_ids,:);
      rand_ns(rand_ids,:)=[];
      %%%search over the ns sites and randomly place a vacany at a nn site --  done in nv section
    end
    
    rand_ns(:,2) = 2; %change type to N
    
    nn = zeros(1,5); %find nn's for ns sites
    n_min = [];
    for i = 1:size(rand_ns,1)
      nn_dist = zeros(size(data,1), 2);
      for j = 1:size(data, 1) %find nn distance to ns site
        nn_dist(j,2) = sqrt((rand_ns(i,3)-data(j,3))^2 + ...
            (rand_ns(i,4)-data(j,4))^2 + (rand_ns(i,5)-data(j,5))^2);
        nn_dist(j,1) = j;
      end
      
      [~,arb] = sort(nn_dist(:,2)); %soft nn list 
      nn_dist = nn_dist(arb,:);
      num_min = sum(nn_dist(:,2)<=CC_len+0.25); %num of 1st nn
      min_index = nn_dist(1:num_min,1); %nn id's
      n_min = [n_min; num_min]; %record the num of 1st nn's. Useful if not all sites have same num 1st nn's

      nn = [nn; data(min_index,:)]; 
      data(min_index,:) = []; %save nn's to shift and delete from "data" so vac.'s wont be placed there later
      
    end
    nn(1,:) = []; %nn data 
    tmp_nn = zeros(1,5); 
    
    %%% now we stretch the bonds along each nn-ns axis by moving the nn sites
    for k = 1:size(rand_ns, 1)
      for l = 1:n_min(k)
        
        if nn(l,3) < rand_ns(k,3)
          %subtract 0.25 from x coordinates - calculated elsewhere. changed
          %cc - bond len to n-c bond len = 1.28%
          nn(l,3) = nn(l,3) - d_nc;
        else  
          %add 0.25 to x coordinates
          nn(l,3) = nn(l,3) + d_nc;
        end
        
        if nn(l,4) < rand_ns(k,4)
          %subtract 0.25 from y coordinates
          nn(l,4) = nn(l,4) - d_nc;
        else  
          %add 0.25 to y coordinates
          nn(l,4) = nn(l,4) + d_nc;
        end
        
        if nn(l,5) < rand_ns(k,5)
          %subtract 0.25 from z coordinates
          nn(l,5) = nn(l,5) - d_nc;
        else  
          %add 0.25 to z coordinates
          nn(l,5) = nn(l,5) + d_nc;
        end
                         
      end  
      tmp_nn = [tmp_nn; nn(1:n_min(k),:)];
      nn(1:n_min(k),:)=[];
      
    end
    tmp_nn(1,:) = [];
    nn_ns = tmp_nn;
  else
   rand_ns = [];
   nn_ns = [];
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if nv ~= 0
   %stretch bonds for v-ns to 1.69 and v-c to 1.63 and then delete the v
   %find nn's for n sites to remove a site for v and stretch bonds
   %Actually, the v is taken to be the point chosen first and then a random
   %nn is set as a ns site. The reason is it is easier to strethc the bonds
   %centered around the v site.
    nv_nn = zeros(1,5); 
    nv_min = [];
    for i = 1:size(rand_nv,1)
      nn_dist = zeros(size(data,1), 2);
      for j = 1:size(data, 1) %find nn distance to v site
        nn_dist(j,2) = sqrt((rand_nv(i,3)-data(j,3))^2 + ...
            (rand_nv(i,4)-data(j,4))^2 + (rand_nv(i,5)-data(j,5))^2);
        nn_dist(j,1) = j;
      end
      
      [~,arb] = sort(nn_dist(:,2)); %sort nn list 
      nn_dist = nn_dist(arb,:);
      num_min = sum(nn_dist(:,2)<=CC_len+0.25); %num of 1st nn
      min_index = nn_dist(1:num_min,1); %nn id's
      nv_min = [nv_min; num_min]; %record the num of 1st nn's. Useful because not all sites have same num 1st nn's
      % i.e. sites on the edge or corner do not have 4 nn's

      nv_nn = [nv_nn; data(min_index,:)]; 
      data(min_index,:) = []; %save nn's to shift and delete from "data" so other vac.'s wont be placed there later
      
    end
    nv_nn(1,:) = []; %nn data 
    tmp_nn = zeros(1,5);
    
    for p = 1:size(rand_nv, 1) %find nn's for site p
        tmp_c = zeros(1,5);
        tmp_c = [tmp_c; nv_nn(1:nv_min(p),:)];
        tmp_c(1,:) = [];
        
        rand_id = randperm(size(tmp_c,1), 1);
        n_site = tmp_c(rand_id,:); %isolate one site as ns, rest are Carbon
        n_site(2) = 2;
        tmp_c(rand_id,:) = []; %remove n site from carbon nn's
        
        for q = 1:size(tmp_c, 1) %loop over c nn sites to stretch bonds
            if tmp_c(q,3) < rand_nv(p,3)
              %subtract 0.0486 from x coordinates - calculated elsewhere. changed
              %cc - bond len to nv - c bond len = 1.28%
              tmp_c(q,3) = tmp_c(q,3) - d_nvc;
            else  
              %add 0.0486 to x coordinates
              tmp_c(q,3) = tmp_c(q,3) + d_nvc;
            end

            if tmp_c(q,4) < rand_nv(p,4)
              %subtract 0.0486 from y coordinates
              tmp_c(q,4) = tmp_c(q,4) - d_nvc;
            else  
              %add 0.0486 to y coordinates
              tmp_c(q,4) = tmp_c(q,4) + d_nvc;
            end

            if tmp_c(q,5) < rand_nv(p,5)
              %subtract 0.0486 from z coordinates
              tmp_c(q,5) = tmp_c(q,5) - d_nvc;
            else  
              %add 0.0486 to z coordinates
              tmp_c(q,5) = tmp_c(q,5) + d_nvc;
            end
        end
        
        if n_site(3) < rand_nv(p,3) %stretch n to v distance
          %subtract 0.0833 from x coordinates - calculated elsewhere. changed
          %cc - bond len to nv - c bond len = 1.28%
          n_site(3) = n_site(3) - d_nvn;
        else  
          %add 0.0833 to x coordinates
          n_site(3) = n_site(3) + d_nvn;
        end

        if n_site(4) < rand_nv(p,4)
          %subtract 0.0833 from y coordinates
          n_site(4) = n_site(4) - d_nvn;
        else  
          %add 0.0833 to y coordinates
          n_site(4) = n_site(4) + d_nvn;
        end

        if n_site(5) < rand_nv(p,5)
          %subtract 0.0833 from z coordinates
          n_site(5) = n_site(5) - d_nvn;
        else  
          %add 0.0833 to z coordinates
          n_site(5) = n_site(5) + d_nvn;
        end
        
        tmp_nn = [tmp_nn; tmp_c; n_site];
    end 
    
    tmp_nn(1,:) = [];
    nv_nn = tmp_nn;
    
  else
      nv_nn = []; %if no nv centers, set nv_nn to empty
  end
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  if v ~= 0
    rand_ids = randperm(size(data, 1), v); %randomly choose some of the ns centers to form nv pairs
    rand_v = data(rand_ids,:);
    data(rand_ids,:) = []; 
    %%% deleting the randomly selected v atoms from "data" is
    %%% equivalent to "placing a vacancy"
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if v2 ~= 0
      
    v2_ids = randperm(size(rand_v,1), v2); %randomly choose v2 sites from v
    v2_sites = rand_v(v2_ids,:);
    rand_v(v2_ids,:) = []; %not necessary to delete v2 sites from v...
    
    nn = zeros(v2,5);
    for i = 1:size(v2_sites,1) %find nn's to v2 sites
      nn_dist = zeros(size(data,1), 2);
      for j = 1:size(data, 1)
        nn_dist(j,2) = sqrt((v2_sites(i,3)-data(j,3))^2 + ...
            (v2_sites(i,4)-data(j,4))^2 + (v2_sites(i,5)-data(j,5))^2);
        nn_dist(j,1) = j;
      end
      
      [~,arb] = sort(nn_dist(:,2));
      nn_dist = nn_dist(arb,:);
      num_min = sum(nn_dist(:,2)<=CC_len+0.25);
      min_id = nn_dist(1:num_min,1); %nn id's
      
      rand_id = randperm(size(min_id, 1), 1); %randomly pick 1 nn as v2 site
      min_id = min_id(rand_id,1);
      nn(i,:) = data(min_id,:); %record v2 data just cause
      data(min_id, :) = []; %remove v2 site from data

    end
    
    v_nn = nn; %just to record ids
    
  end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    Pos = [data; rand_ns; nn_ns; nv_nn]; %add ns, shifted ns-nn, and shifter nv-nn sites back to data
    [~,arb] = sort(Pos(:,3));
    Pos = Pos(arb, :); %sort Pos by x-coord
    
    save('tmp.mat');
   
    if lammps == 1
        LAMMPS(name, Pos);
    end

    if xyz == 1
        XYZ(name, Pos);
    end
end


function LAMMPS(name, Pos)
    buffer = 0.009;
    masses = [12.0107; 14.0067];
    xlo = min(Pos(:,3))-buffer; xhi = max(Pos(:,3))+buffer;
    ylo = min(Pos(:,4))-buffer; yhi = max(Pos(:,4))+buffer;
    zlo = min(Pos(:,5))-buffer; zhi = max(Pos(:,5))+buffer;
    sys = size(Pos, 1);

    %% Build the Header
    fid = fopen(sprintf('data.NV.%s',name),'w');
    fprintf(fid,'LAMMPS data file\r\n\r\n%d atoms\r\n\r\n%d atom types\r\n\r\n'...
        ,sys,length(masses));
    fprintf(fid,'%6.3f %6.3f xlo xhi\r\n%6.3f %6.3f ylo yhi\r\n%6.3f %6.3f zlo zhi\r\n'...
        ,xlo,xhi,ylo,yhi,zlo,zhi);
    fprintf(fid,'\r\nMasses\r\n\r\n');
    for i = 1:length(masses)
        fprintf(fid,'%d %6.4f\r\n',i,masses(i));
    end
    fprintf(fid,'\r\nAtoms\r\n\r\n');
    %% Atoms
    for i = 1:sys  
        fprintf(fid,'%d\t%d\t%6.3f\t%6.3f\t%6.3f\r\n',i,Pos(i,2),Pos(i,3),Pos(i,4),Pos(i,5));
    end
    fclose(fid);  
end


function XYZ(name ,Pos)
    [~,arb] = sort(Pos(:,2)); %sort by atom type
    Pos = Pos(arb,:);
    sys = size(Pos, 1);
    fid = fopen(sprintf('NV.%s.xyz',name),'w');
    fprintf(fid,'%d\r\n',sys);
    fprintf(fid,'%4.2f\t0\t0\t0\t%4.2f\t0\t0\t0\t%4.2f\r\n',max(Pos(:,3)),max(Pos(:,4))...
        ,max(Pos(:,5)));
    for i = 1:sys
        if Pos(i,2) == 1
            fprintf(fid,'C\t%4.2f\t%4.2f\t%4.2f\r\n',Pos(i,3),Pos(i,4),Pos(i,5));
        elseif Pos(i,2) == 2
            fprintf(fid,'N\t%4.2f\t%4.2f\t%4.2f\r\n',Pos(i,3),Pos(i,4),Pos(i,5));
        end
    end
    fclose(fid);
end