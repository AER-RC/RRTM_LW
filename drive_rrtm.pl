#!/usr/bin/perl -w


# This driver runs various versions of RRTM for a series of cases and stores the output

$band_gen_dir="/project/p1905/band_build_with_flux/test_gb3/";

$code = ("/project/p1905/band_build/rrtm_lw_so2_kaust/rrtm_v3.3.1_linux_pgf90");

$band_id = '3';
$band_num = $band_id;
$band_name = 'band'.$band_id;

$rrtm_dir = './test_'.$band_name;

$kg_file = 'k_gB0'.$band_id.'.f';
$kg_file_base = $kg_file.'.base';
print $kg_file,$kg_file_base,"\n";

print $band_num,"\n";

chdir($rrtm_dir);
@case_names =  glob("INPUT_RRTM.GARAND*o*");
$ncases = @case_names;
print $ncases,"\n";

@case_tags = ();
for ($ic=0;$ic<$ncases; $ic++)  {
   $case_tags[$ic] = substr($case_names[$ic],11);
}
chdir("../");

@kminor_files = glob("$band_gen_dir/k_minor*.txt");
$nfiles = @kminor_files;
print $nfiles,"\n";

foreach $il (@kminor_files) {
   chdir ("src/");
   system ("/bin/rm $kg_file");
   system ("/bin/cp $kg_file_base $kg_file");
   system ("cat $il >> $kg_file");

   $minor_map = substr(`basename $il`,8);
   $minor_map = substr($minor_map,0,-5);
   print $minor_map,"\n";
   system ("cp $kg_file $kg_file.$minor_map");

   chdir ("../");
   system ("make -f makefiles/make_rrtm");
   system ("ls -l rrtm_v3.3.1_linux_pgf90");

   chdir($rrtm_dir);
   for ($ic=0;$ic<$ncases;$ic++)  {
      system ("/bin/rm INPUT_RRTM OUTPUT_RRTM");
      system ("cp $case_names[$ic] INPUT_RRTM");
      system ($code);
      print ($case_tags[$ic],"\n");
      print ($band_name,"\n");
      print ($minor_map,"\n");
      system ("/bin/mv OUTPUT_RRTM OUTPUT_RRTM.$case_tags[$ic].$band_name.$minor_map");
   }
   chdir ("../");
}
