pro gen_runs, freq, submit=submit, ifield_start=ifield_start, ifield_end=ifield_end, partition=partition

    num_fields = 20
    num_fields_node = 2
    num_cpus_field  = 8

    freq_str = strcompress(string(freq),/remove)
    ;dr_str = strcompress(string(planck_dr),/remove)
    dr_str= 'R2.00'
    
    f = lps12_fieldstruct()
    
    if (keyword_set(ifield_start)) then ifield0=ifield_start else ifield0=0
    if (keyword_set(ifield_end)) then ifield1=ifield_end else ifield1=num_fields-1
    if (keyword_set(partition)) then part=partition else part='kicp'

    for ifield=ifield0, ifield1 do begin
        field_str = strcompress(string(ifield),/remove)
        field_name = f[ifield].name
        submit_file = 'scripts/submit_hfi_'+freq_str+'_'+dr_str+'_'+field_name+'.sh'

        openw, 5, submit_file
        printf, 5, '#!/bin/bash'
        printf, 5, '#SBATCH --partition='+part
        printf, 5, '#SBATCH --account=kicp'
        printf, 5, '#SBATCH --nodes=1'
        printf, 5, '#SBATCH --ntasks=1'
        printf, 5, '#SBATCH --cpus-per-task='+strcompress(string(num_cpus_field),/remove)
        printf, 5, '#SBATCH --job-name='+field_name+'_'+freq_str
        printf, 5, '#SBATCH --time=36:00:00'
        printf, 5, ' '
        printf, 5, 'export HOSTNAME=midway'
        printf, 5, 'export DATAPATH=/home/zhenhou/data/'
        ;printf, 5, 'source /home/zhenhou/software/idl/idl/bin/idl_setup.bash'
        printf, 5, 'export ITT_DIR=/home/zhenhou/software/idl'
        printf, 5, 'export IDL_DIR=/home/zhenhou/software/idl/idl70'
        printf, 5, 'cd /home/zhenhou/Projects/projects/spt_x_planck/planck_sims/powspec/'
        printf, 5, ' '
        printf, 5, 'export IDL_PATH IDL_STARTUP'
        printf, 5, 'OIDL_PATH="${IDL_PATH-<IDL_DEFAULT>}"'
        printf, 5, 'OIDL_STARTUP="${IDL_STARTUP}"'
        printf, 5, 'IDL_PATH="+/home/zhenhou/Projects/spt_cvs/spt_analysis:+/home/zhenhou/Projects/CMBtools/idl_routines:+${HEALPIX}/src/idl:${OIDL_PATH}"'
        printf, 5, 'IDL_STARTUP="${HEALPIX}/src/idl/HEALPix_startup"'
        printf, 5, 'echo $IDL_STARTUP'
        printf, 5, '$IDL_DIR/bin/idl -e "make_half_sims_xspec, '+field_str+', '+freq_str+'"'
        close, 5

        if (keyword_set(submit)) then begin
            spawn, 'sbatch '+submit_file
        endif
    endfor

end
