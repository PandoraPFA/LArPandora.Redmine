#!/usr/bin/perl

# This program creates GDML sub-files, with values supplied by user
# parameters.  Geometry/gdml/make_gdml.pl "zips" together those
# sub-files to make a single detector description.

# Packages
use Math::Trig;
use XML::LibXML;
use Getopt::Long;

# Get the input parameters from an XML file. Optionally append a
# suffix to the GDML sub-files we create.

GetOptions( "input|i:s" => \$input,
	    "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output);

if ( defined $help )
{
    # If the user requested help, print the usage notes and exit.
    usage();
    exit;
}

if ( ! defined $suffix )
{
    # The user didn't supply a suffix, so append nothing to the file
    # names.
    $suffix = "";
}
else
{
    # Otherwise, stick a "-" before the suffix, so that a suffix of
    # "test" applied to filename.gdml becomes "filename-test.gdml".
    $suffix = "-" . $suffix;
}

# Create an XML parser.
$parser = new XML::LibXML;

# Read in the parameters from an XML file. The following command
# slurps the entire file into a DOM data structure.
$xmldata = $parser->parse_file($input);

# Go through each parameter in the DOM data structure:
foreach $parameter ( $xmldata->findnodes('/parameters/geometry/parameter') )
{
    # Get the name and value attributes for that parameter:
    $name = $parameter->getAttribute("name");
    $value = $parameter->getAttribute("value");

    # Here's the clever part: The following eval creates a variable
    # with the same name as $name. For example, if $name eq "TPCDepth",
    # then the following statement assigns the value to $TPCDepth. The
    # value is in quotes, because some of the parameters have text
    # strings in them (like "kInch").

    eval "\$$name = '$value'";
}

# Our calculations and constants depend on the geometry of the wires.
$SinUVAngle = sin( deg2rad($UVAngle) );
$CosUVAngle = cos( deg2rad($UVAngle) );
$TanUVAngle = tan( deg2rad($UVAngle) );

# The routines that create the GDML sub-files. Most of the explanatory
# comments are in gen_defs().
$wires_on=1; 			# turn wires on=1 or off=0
$wire_int=20;
$NumberOfTPCPlanes=3;
$tpc_neg_length=1000;
$tpc_neg_width=200;
$tpc_neg_height=200;
$wire_frame_width=9.5;
$wire_plane_height=240;
$wire_plane_length=1042;
$wires_plength=($wire_plane_length - 2*$wire_frame_width) ;
$wires_pwidth=($wire_plane_height - 2*$wire_frame_width) ;
$field_cage_width		=	200;
$field_cage_height		=	180;
$field_cage_cross_length	=	sqrt(($field_cage_width)**2+($field_cage_height-50)**2);
$field_cage_length		=	1000;
$field_cage_loop_interval	=	1; 	# =1 is normal, =4 skips 3/4
$spacers_on_off		= 	"off"; 	# "on" or "off" for tube spacers (off saves time)
$electronics_height=12;
$pmt_switch="off";		#turn on or off depending on pmts wanted


gen_defs();
gen_rotations();
gen_materials();

gen_tpcplanevert();
gen_tpcplane();

 gen_groundplate();	# physical volumes defined in gen_tpc()
 gen_cathode();		# physical volumes defined in gen_tpc()
 gen_fieldcage();	# physical volumes defined in gen_tpc()
gen_tpc();

if ( $pmt_switch eq "on" ) {  gen_pmt();	}	# physical volumes defined in gen_cryostat()
gen_cryostat();

gen_enclosure();
gen_world();
write_fragments();

exit;



sub usage()
{
    print "Usage: $0 [-h|--help] -i|--input <parameters-file> [-o|--output <fragments-file>] [-s|--suffix <string>]\n";
    print "       -i/--input can be omitted; <parameters-file> contains geometry and material parameters\n";
    print "       if -o is omitted, output goes to STDOUT; <fragments-file> is input to make_gdml.pl\n";
    print "       -s <string> appends the string to the file names; useful for multiple detector versions\n";
    print "       -h prints this message, then quits\n";
}



# Create the detector constant file. This file is actually temporary,
# since the make_gdml.pl program will interpret its contents rather
# than include it in the final GDML file.

sub gen_defs()
{
    # Set up the output file.
    $CONSTANTS = "microboone/micro-defs" . $suffix . ".gdml";
    push (@gdmlFiles, $CONSTANTS); # Add file to list of constant files
    $CONSTANTS = ">" . $CONSTANTS;
    open(CONSTANTS) or die("Could not open file $CONSTANTS for writing");

    # Create some math constants.
    my $pi = pi;

    # Though it's not strictly necessary, make each sub-file a valid
    # XML (if not GDML) document; that makes it accessible to an XML
    # parser if needed.

    # Here is a neat way to print out a block of text without getting
    # involved in a lot of messy quoting and formatting with print
    # statements.

    print CONSTANTS <<EOF;
<?xml version='1.0'?>
<define>
<constant name="kInch"	value="2.54" />
<constant name="kPi"	value="$pi" />
<constant name="kDetEnclosureWidth"	  value="$DetEnclosureWidth" />
<constant name="kDetEnclosureHeight"	  value="$DetEnclosureHeight" />
<constant name="kDetEnclosureLength"	  value="$DetEnclosureLength" />
<constant name="kDirtThickness"           value="$DirtThickness" />
<constant name="kWorldW"                  value="100.0*kDetEnclosureWidth"/>
<constant name="kWorldH"                  value="100.0*kDetEnclosureHeight"/>
<constant name="kWorldL"                  value="100.0*kDetEnclosureLength"/>

<constant name="kTPCWidth"		  value="$TPCWidth" />
<constant name="kTPCLength"		  value="$TPCLength" />
<constant name="kTPCDepth"		  value="$TPCDepth" />
<constant name="kTPCWallThickness"        value="$TPCWallThickness" />

<constant name="kTPCWirePlaneThickness"   value="$TPCWirePlaneThickness" />
<constant name="kTPCWireThickness"	  value="$TPCWireThickness" />
<constant name="kTPCWirePlaneWidth"       value="$wires_pwidth" />
<constant name="kTPCWirePlaneLength"	  value="$wires_plength" />
<constant name="kWireFrameDepth"  	  value="9" /> 
<constant name="kWireFrameWidth"  	  value="$wire_frame_width" /> 
<constant name="kWirePlaneHeight"  	  value="$wire_plane_height" /> 
<constant name="kWirePlaneLength"  	  value="$wire_plane_length" /> 
<constant name="kWireFrameVInHeight"  	  value="0.5*(kWirePlaneHeight-3*kWireFrameWidth)" /> 

<constant name="kTPCWirePitch"            value="$TPCWirePitch"/>
<constant name="kSinUVAngle"              value="$SinUVAngle"/>
<constant name="kCosUVAngle"              value="$CosUVAngle"/>
<constant name="kTanUVAngle"              value="$TanUVAngle"/>
<constant name="kTPCWireXPitch"           value="kTPCWirePitch/kCosUVAngle"/>

<constant name="kCathodeFrameWidth"  	value="9" /> 
<constant name="kCathodeFrameDepth"  	value="5" /> 
<constant name="kCathodePlateDepth"  	value="0.1" /> 
<constant name="kCathodeWidth"  	value="5.1" /> 
<constant name="kCathodeHeight"  	value="240" /> 
<constant name="kCathodeLength"  	value="1042" /> 
<constant name="kCathodeFrameVInHeight"  	value="0.5*(kCathodeHeight-3*kCathodeFrameWidth)" /> 

<constant name="kGroundPlateWidth"      value="224" /> 
<constant name="kGroundPlateHeight"     value="0.1" /> 
<constant name="kGroundPlateLength"             value="1100" /> 
<constant name="kGroundBeamWidth"       value="2.5" /> 
<constant name="kGroundBeamHeight"      value="2.5" /> 
<constant name="kGroundBeamLength"              value="kGroundPlateLength" />
<constant name="kGroundBeamThickness"   value=".15" /> 
</define>
EOF

   close(CONSTANTS);
}


sub gen_rotations()
{
    my $WirePlusRotation = $UVAngle + 90;
    my $WireMinusRotation = $UVAngle - 90;

    $ROTATIONS = "microboone/micro-rotations" . $suffix . ".gdml";
    push (@gdmlFiles, $ROTATIONS); # Add file to list of GDML fragments
    $ROTATIONS = ">" . $ROTATIONS;
    open(ROTATIONS) or die("Could not open file $ROTATIONS for writing");

    print ROTATIONS <<EOF;
<?xml version='1.0'?>
<define>
   <rotation name="rPlus30AboutX"  unit="deg" x="30"  y="0"   z="0"/>
   <rotation name="rPlus60AboutX"  unit="deg" x="60"  y="0"   z="0"/>
   <rotation name="rPlus90AboutX"  unit="deg" x="90"  y="0"   z="0"/>
   <rotation name="rMinus90AboutX"  unit="deg" x="-90"  y="0"   z="0"/>
   <rotation name="rPlusUVAngleAboutX"  unit="deg" x="150" y="0"   z="0"/>
   <rotation name="rPlus150AboutX"      unit="deg" x="150" y="0"   z="0"/>
   <rotation name="rPlus180AboutX"      unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rMinusUVAngleAboutX" unit="deg" x="-30" y="0"   z="0"/>
   <rotation name="rPlus30AboutY"  unit="deg" x="0"   y="30"  z="0"/>
   <rotation name="rPlus60AboutY"  unit="deg" x="0"   y="60"  z="0"/>
   <rotation name="rPlus90AboutY"  unit="deg" x="0"   y="90"  z="0"/>
   <rotation name="rPlus180AboutY" unit="deg" x="0"   y="180" z="0"/>
   <rotation name="rMinus90AboutY" unit="deg" x="0"   y="-90" z="0"/>
   <rotation name="rPlus90AboutZ"  unit="deg" x="0"   y="0"   z="90"/>
   <rotation name="rMinus90AboutZ"  unit="deg" x="0"   y="0"   z="-90"/>
   <rotation name="rPlus180AboutZ"      unit="deg" x="0"   y="0"   z="180"/>
   <rotation name="rMinus180AboutZ"     unit="deg" x="0"   y="0"   z="-180"/>
   <rotation name="rMinus90AboutYPlus180AboutZ" unit="deg" x="0" y="-90" z="180"/>
   <rotation name="rMinus90AboutYMinus90AboutZ" unit="deg" x="0" y="-90" z="-90"/>
   <rotation name="rPlus90AboutYPlus180AboutZ" unit="deg" x="0" y="90" z="180"/>
   <rotation name="rMinus90AboutYPlus90AboutZ" unit="deg" x="0" y="-90" z="90"/>
   <rotation name="rPlus90AboutYMinus90AboutZ" unit="deg" x="0" y="90" z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutZ"  unit="deg" x="90" y="0"   z="90"/>
   <rotation name="rPlus90AboutXPlus180AboutZ" unit="deg" x="90" y="0"   z="180"/>
   <rotation name="rPlus90AboutXMinus90AboutY" unit="deg" x="90" y="-90" z="0"/>
   <rotation name="rPlus90AboutXMinus90AboutZ" unit="deg" x="90" y="0"   z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutY"  unit="deg"  x="90" y="90" z="0"/>
   <rotation name="rPMTRotation1"  unit="deg" x="90"  y="270"   z="0"/>
</define>
EOF
    close (ROTATIONS);
}


sub gen_materials()
{
    # Create the materials file name and open it.
    $MATERIALS = "materials" . $suffix . ".gdml";
    push (@gdmlFiles, $MATERIALS); # Add file to list of GDML fragments
    $MATERIALS = ">" . $MATERIALS;
    open(MATERIALS) or die("Could not open file $MATERIALS for writing");

    # Write the standard XML prefix.
    print MATERIALS <<EOF;
<?xml version='1.0'?>
EOF

    # Go back the DOM structure read in near the beginning of the
    # program. For each <materials /> element (and there'll probably
    # be only one):
    foreach $materials ( $xmldata->findnodes('/parameters/materials') )
    {
	# Convert that element back to text, and write it out.
	print MATERIALS $materials->toString;
    }

    close (MATERIALS);
}


# This is a re-write of Brian Rebel's gen_microvertplane.C into
# Perl. It contructs the TPC wire plane for the Y view.

sub gen_tpcplanevert()
{


##### temporary edits:
#####   - TPCPlaneVert   y="TPCWidth" z="kTPCLength"
#####   - TPCWireVert    y="TPCWidth"
#####   - (above)my $NumberWires = int($TPCLength / $TPCWirePitch ) - 1



    my $NumberWires = int( ( $wires_plength ) / $TPCWirePitch ) - 1;

    $GDML = "microboone/micro-vertplane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Define the solids and structures: the wires and the TPC wire plane.

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
<tube name="TPCWireVert"
  rmax="0.5*kTPCWireThickness"
  z="$wires_pwidth"		
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
<box name="TPCPlaneVert"
  x="kTPCWirePlaneThickness" 
  y="$wires_pwidth" 
  z="$wires_plength"
  lunit="cm"/>
</solids>
<structure>
  <volume name="volTPCWireVert">
    <materialref ref="Titanium"/>
    <solidref ref="TPCWireVert"/>
  </volume>
  <volume name="volTPCPlaneVert">
    <materialref ref="LAr"/>       
    <solidref ref="TPCPlaneVert"/>
EOF

    # the wires 
    for ( $i = 0; $i < ( $NumberWires / $wire_int ); ++$i )
    {
	$j=($wire_int*$i);
	print GDML <<EOF;
	    <physvol>
	     <volumeref ref="volTPCWireVert"/>
	     <position name="posTPCWireVert$i" unit="cm" z="(-0.5*kTPCWirePlaneLength)+kTPCWirePitch*($j+1)" x="0" y="0"/>
	     <rotationref ref="rPlus90AboutX"/>
	    </physvol>
EOF
    }

    print GDML <<EOF;
  </volume>
</structure>
</gdml>
EOF

    close(GDML);
}

# This is a re-write of Brian Rebel's gen_microplane.C into Perl. It
# constructs the TPC wire plane for the U or V view.

sub gen_tpcplane()
{

#### temporary edits
####    - my $NumberWires = $TPCLength / $TPCWirePitch - 1;
####    - my $NumberWiresPerEdge = int( $TPCLength / $TPCYWirePitch );
####    - my $NumberSideWires = int( $TanUVAngle * $TPCWidth / $TPCYWirePitch );
####    - <tube name="TPCWireCommon" rmax="0.5*kTPCWireThickness" z="kTPCWidth/kCosUVAngle" deltaphi="2*kPi" aunit="rad" lunit="cm"/>
####    - <box name="TPCPlane"  x="kTPCWirePlaneThickness"  y="kTPCWidth"  z="kTPCLength"  lunit="cm"/>

    my $NumberWires = ( $wires_plength ) / $TPCWirePitch - 1;

    $GDML = "microboone/micro-plane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Calculate the number of wire ends on a given y-edge of the plane.
    my $TPCYWirePitch = $TPCWirePitch / $CosUVAngle;
    my $NumberWiresPerEdge = int( ( $wires_plength ) / $TPCYWirePitch );

    # How many side wires will be "cut off" by the lower or higher
    # z-edge?
    my $NumberSideWires = int( $TanUVAngle * ( $wires_pwidth ) / $TPCYWirePitch );

    # The number of full-length "center" wires.
    my $NumberCenterWires = $NumberWiresPerEdge - $NumberSideWires;

    # define the solids
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
EOF

    # wires on either end of the tpc
    for($i = 0; $i < $NumberSideWires; ++$i)
    {
	print GDML <<EOF;
<tube name="TPCWire$i"
  rmax="0.5*kTPCWireThickness"
  z="kTPCWireXPitch*($i+1)/kSinUVAngle"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
EOF
    }

    # The solids for the middle wire and the TPC wire plane, and start off the structures.
    print GDML <<EOF;
<tube name="TPCWireCommon"
  rmax="0.5*kTPCWireThickness"
  z="($wires_pwidth)/kCosUVAngle"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
<box name="TPCPlane"
  x="kTPCWirePlaneThickness"
  y="($wires_pwidth)"
  z="($wires_plength)"
  lunit="cm"/>
</solids>
<structure>
EOF
 
    # the wires at either end of the plane
    for ($i = 0; $i < $NumberSideWires; ++$i)
    {
	print GDML <<EOF;
    <volume name="volTPCWire$i">
	<materialref ref="Titanium"/>
	<solidref ref="TPCWire$i"/>
    </volume>
EOF
    }

  
    # The wires in the middle of the plane, and the plane itself.
    print GDML <<EOF;
    <volume name="volTPCWireCommon">
      <materialref ref="Titanium"/>
      <solidref ref="TPCWireCommon"/>
    </volume>
    <volume name="volTPCPlane">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

  # the wires at the -z end
  for ($i = 0; $i < $NumberSideWires; ++$i)
  {
      print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWire$i"/> 
     <position name="posTPCWire$i" unit="cm" y="(-0.5*kTPCWirePlaneWidth)+0.5*($i+1)*kTPCWireXPitch/kTanUVAngle" z="-0.5*kTPCWirePlaneLength+0.5*kTPCWireXPitch*($i+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  }

  # The wires in the middle.
  for ($i = 0; $i < $NumberCenterWires - 1 ; ++$i)
  {
      my $j = $NumberSideWires+$i;
      print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWireCommon"/>
     <position name="posTPCWire$j" unit="cm" y="0" z="(-0.5*kTPCWirePlaneLength)+kTPCWireXPitch*(0.5*$NumberSideWires + $i+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  }

  # the wires at the +z end
  for ($i = 0; $i < $NumberSideWires; ++$i)
  {
      my $j = $NumberSideWires-$i-1;
      my $k = $NumberCenterWires+$NumberSideWires+$i;

      print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWire$j"/>
     <position name="posTPCWire$k" unit="cm" y="0.5*kTPCWirePlaneWidth-0.5*($j+1)*kTPCWireXPitch/kTanUVAngle" z="0.5*kTPCWirePlaneLength-0.5*kTPCWireXPitch*($j+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  }

      print GDML <<EOF;
  </volume>
</structure>
</gdml>
EOF

  close(GDML);
}


#
# subdirectory to write field cage
sub gen_fieldcage() {

    # Set up the output file.
    $FIELDCAGE = "microboone/micro-fieldcage.gdml";
    push (@gdmlFiles, $FIELDCAGE); # Add file to list of constant files
    $FIELDCAGE = ">" . $FIELDCAGE;
    open(FIELDCAGE) or die("Could not open file $FIELDCAGE for writing");

    # Print the Field Cage constants
    print FIELDCAGE <<EOF;
<define>
 <constant name="kFieldCageTPCClearance"    	value="5*kInch" />

 <constant name="kFieldCageTubeRadius"  	value="0.5*kInch" />
 <constant name="kFieldCageTubeThickness" 	value="0.25*kInch" />
 <constant name="kFieldCageBeamDepth"           value="12.5"/>
 <constant name="kFieldCageBeamWidth"           value="2.5"/>
 <constant name="kFieldCageCrossDepth"          value="2.5"/>
 <constant name="kFieldCageCrossWidth"          value="4"/>
 <constant name="kFieldCageCrossLength"         value="$field_cage_cross_length"/>

 <constant name="kTPCTotalLength"         	value="$field_cage_length"/>
 <constant name="kTPCTotalWidth"         	value="$field_cage_width"/>
 <constant name="kTPCTotalHeight"         	value="$field_cage_height"/>

 <constant name="kFieldCageLoopLength"       	value="kTPCTotalLength+2*(kFieldCageTPCClearance+2*kFieldCageTubeRadius)"/>
 <constant name="kFieldCageLoopWidth"   	value="$field_cage_width"/>
 <constant name="kFieldCageLoopHeight"       	value="kTPCTotalHeight+2*(kFieldCageTPCClearance+2*kFieldCageTubeRadius)"/>

 <constant name="kFieldCageCornerRadius"	value="0.5*kFieldCageTPCClearance"/>
 <constant name="kFieldCageCornerY"       	value="(0.5*kFieldCageLoopHeight)-kFieldCageCornerRadius-kFieldCageTubeRadius"/>
 <constant name="kFieldCageCornerZ"       	value="(0.5*kFieldCageLoopLength)-kFieldCageCornerRadius-kFieldCageTubeRadius"/>

 <constant name="kFieldCageHeight"       	value="kFieldCageLoopHeight+2*(kFieldCageBeamDepth+kFieldCageCrossDepth)"/>
 <constant name="kFieldCageLength"       	value="kFieldCageLoopLength+2*(kFieldCageBeamDepth+kFieldCageCrossDepth)"/>
 <constant name="kFieldCageWidth"      		value="$field_cage_width"/>

 <constant name="kFieldCageBeamYInt"            value="0.5*(kFieldCageLoopHeight-50)"/>
 <constant name="kFieldCageBeamZPos"            value="0.5*(kFieldCageLoopLength)"/>
 <constant name="kFieldCageBeamYPos"            value="0.5*(kFieldCageLoopHeight)"/>
 <constant name="kFieldCageBeamZInt"            value="0.5*(kFieldCageLoopLength-40)"/>

 <constant name="kFieldCageCrossYPos"           value="0.5*(kFieldCageLoopHeight+kFieldCageCrossDepth)+kFieldCageBeamDepth"/>
 <constant name="kFieldCageCrossZPos"           value="0.5*(kFieldCageLoopLength+kFieldCageCrossDepth)+kFieldCageBeamDepth"/>
</define>
EOF

    # Prints Field Cage solids
    print FIELDCAGE <<EOF;
<solids>
 <tube name="FieldCageTubeZ"
  rmin="kFieldCageTubeRadius - kFieldCageTubeThickness"
  rmax="kFieldCageTubeRadius"  
  z="kFieldCageLoopLength-2*(kFieldCageCornerRadius+kFieldCageTubeRadius)" 
  deltaphi="2*kPi" 
  aunit="rad" 
  lunit="cm"/> 
 <tube name="FieldCageTubeY" 
  rmin="kFieldCageTubeRadius - kFieldCageTubeThickness"
  rmax="kFieldCageTubeRadius"  
  z="kFieldCageLoopHeight-2*(kFieldCageCornerRadius+kFieldCageTubeRadius)"  
  deltaphi="2*kPi" 
  aunit="rad" 
  lunit="cm"/> 
</solids> 
EOF

    # Prints Field Cage tube loop sub-structure
    print FIELDCAGE <<EOF;
<structure> 
 <volume name="volFieldCageTubeTop"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeZ"/> 
 </volume> 
 <volume name="volFieldCageTubeBot"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeZ"/> 
 </volume> 
 <volume name="volFieldCageTubeFront"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeY"/> 
 </volume> 
 <volume name="volFieldCageTubeBack"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeY"/> 
 </volume> 
</structure>
EOF

close(FIELDCAGE); 

}


# Cathode solids and volumes
sub gen_cathode() {

    $CATHODE = "microboone/micro-cathode" . $suffix . ".gdml";
    push (@gdmlFiles, $CATHODE); # Add file to list of GDML fragments
    $CATHODE = ">" . $CATHODE;
    open(CATHODE) or die("Could not open file $CATHODE for writing");

    print CATHODE <<EOF;

<?xml version='1.0'?>
<gdml>
<solids>
 <box name="CathodePlate"
  lunit="cm"
  x="kCathodePlateDepth"
  y="kCathodeHeight"
  z="kCathodeLength"/>
</solids>
<structure>
 <volume name="volCathodePlate">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodePlate"/>
 </volume>
</structure>
</gdml>
EOF

}

sub gen_groundplate() {
#this subroutine will produce the gdml fragment for ground plate

    $GROUNDPLATE = "microboone/micro-groundplate" . $suffix . ".gdml";
    push (@gdmlFiles, $GROUNDPLATE); # Add file to list of GDML fragments
    $GROUNDPLATE = ">" . $GROUNDPLATE;
    open(GROUNDPLATE) or die("Could not open file $GROUNDPLATE for writing");

    print GROUNDPLATE <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="GroundPlate"
  lunit="cm"
  x="kGroundPlateWidth"
  y="kGroundPlateHeight"
  z="kGroundPlateLength"/>
</solids>
<structure>
 <volume name="volGroundPlate">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="GroundPlate"/>
 </volume>
</structure>
</gdml>
EOF

    close(GROUNDPLATE);
}


# Parameterize the TPC and the planes within it.
sub gen_tpc()
{
    # Set up the output file.
    $GDML = "microboone/micro-tpc" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="TPC"
  lunit="cm"
  x="kTPCDepth+kTPCWallThickness+(2*kTPCWirePlaneThickness)"
  y="kTPCWidth+(2*kTPCWallThickness)+($electronics_height)"
  z="kTPCLength+(2*kTPCWallThickness)"/>
</solids>
<structure>
 <volume name="volTPC">
  <materialref ref="LAr"/>
 <solidref ref="TPC"/>
EOF

    # Ground Plane physical volumes
    # Center = (9,124,0)
    $ground_plate_X=9;
    $ground_plate_Y=124;
    $ground_plate_Z=0;

    print GDML <<EOF;
  <physvol>
   <volumeref ref="volGroundPlate"/>
   <position name="posGroundPlate" unit="cm" x="$ground_plate_X+0.25" y="$ground_plate_Y-0.5*(kGroundBeamHeight)" z="0"/>
  </physvol>
EOF


    # Cathode Plane physical volumes
    # Center = (0.5*(kTPCWidth-kCathodeWidth),0,0)
    $Cathode_X="0.5*(kTPCWidth-kCathodeWidth)";
    $Cathode_plate_offset="0.5*kCathodePlateDepth";
    $Cathode_frame_position="$Cathode_X+($Cathode_plate_offset)";

    print GDML <<EOF;
  <physvol>
   <volumeref ref="volCathodePlate"/>
   <position name="posCathodePlate" unit="cm" x="$Cathode_X-0.5*(kCathodeFrameDepth)-0.1" y="0" z="0"/>
  </physvol>
EOF


    # Wire Plane physical volumes 
    # Center = ( -0.5*(kTPCWidth-kWireFrameWidth) , 0 , 0 )
    print GDML <<EOF;
  <physvol>
   <volumeref ref="volTPCPlaneVert"/>
   <position name="posTPCPlaneVert" unit="cm" x="-0.5*(kTPCWidth-kWireFrameWidth)-(2*$TPCWirePlaneSpacing)" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCPlane"/>
   <position name="posTPCPlane" unit="cm" x="-0.5*(kTPCWidth-kWireFrameWidth)-$TPCWirePlaneSpacing" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCPlane"/>
   <position name="posTPCPlane2" unit="cm" x="-0.5*(kTPCWidth-kWireFrameWidth)" y="0" z="0"/>
   <rotationref ref="rPlus180AboutY"/>
  </physvol>
EOF


$space=0;
$i=1;
while ( $space < ( $field_cage_width / 2 ) ) {
	$xPos=$space+2;
	print GDML <<EOF;
  <physvol>
   <volumeref ref="volFieldCageTubeTop"/>
   <position name="posFieldCageTubeTopA$i" unit="cm" x="$xPos" y="0.5*(kFieldCageLoopHeight-2*kFieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBot"/>
   <position name="posFieldCageTubeBotA$i" unit="cm" x="$xPos" y="-0.5*(kFieldCageLoopHeight-2*kFieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeFront"/>
   <position name="posFieldCageTubeFrontA$i" unit="cm" x="$xPos" y="0" z="0.5*(kFieldCageLoopLength-2*kFieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertPlusA$i" unit="deg" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBack"/>
   <position name="posFieldCageTubeBackA$i" unit="cm" x="$xPos" y="0" z="-0.5*(kFieldCageLoopLength-2*kFieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertMinusA$i" unit="deg" x="-90" y="0" z="0"/>
  </physvol>
EOF
	print GDML <<EOF;
  <physvol>
   <volumeref ref="volFieldCageTubeTop"/>
   <position name="posFieldCageTubeTopB$i" unit="cm" x="-$xPos" y="0.5*(kFieldCageLoopHeight-2*kFieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBot"/>
   <position name="posFieldCageTubeBotB$i" unit="cm" x="-$xPos" y="-0.5*(kFieldCageLoopHeight-2*kFieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeFront"/>
   <position name="posFieldCageTubeFrontB$i" unit="cm" x="-$xPos" y="0" z="0.5*(kFieldCageLoopLength-2*kFieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertFrontB$i" unit="deg" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBack"/>
   <position name="posFieldCageTubeBackB$i" unit="cm" x="-$xPos" y="0" z="-0.5*(kFieldCageLoopLength-2*kFieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertBackB$i" unit="deg" x="-90" y="0" z="0"/>
  </physvol>
EOF
	$space+=4*$field_cage_loop_interval;
	$i++;
}



    # Closes TPC volume definition space
    print GDML <<EOF;
 </volume> 
</structure>
</gdml>
EOF

   close(GDML);
}


# Generates Ben Jones's PMT micro-pmtdef (with temporary edit to ellipsoid shapes
sub gen_pmt {

    $PMT = "microboone/micro-pmtdef" . $suffix . ".gdml";
    push (@gdmlFiles, $PMT); # Add file to list of GDML fragments
    $PMT = ">" . $PMT;
    open(PMT) or die("Could not open file $PMT for writing");

	print PMT <<EOF;
<solids>
 <tube name="PMTVolume"
  rmax="(6*2.54)"
  z="(11.0*2.54)"
  deltaphi="2*(3.1415926535897)"
  aunit="rad"
  lunit="cm"/>
 <tube name="PMT_TPBCoating"
  rmax="(6.0*2.54)"
  z="0.01"
  deltaphi="2*(3.1415926535897)"
  aunit="rad"
  lunit="cm"/>
 <tube name="PMT_AcrylicPlate"
  rmax="(6.0*2.54)"
  z="(0.2)"
  deltaphi="2*(3.1415926535897)"
  aunit="rad"
  lunit="cm"/>
 <tube name="PMT_Stalk"
  rmax="(1.25*2.54)"
  z="(3.0*2.54)"
  deltaphi="2*(3.1415926535897)"
  aunit="rad"
  lunit="cm"/>
 <tube name="PMT_SteelBase"
  rmax="(6.0*2.54)"
  z="(1.5*2.54)"
  deltaphi="2*(3.1415926535897)"
  aunit="rad"
  lunit="cm"/>
 <tube name="PMT_Underside"
  rmax="2.54*4.0"
  z="2.54*2.5"
  deltaphi="2*3.1415926535897"
  aunit="rad"
  lunit="cm"/>
 <tube name="PMT_Lens"
  rmax="2.54*4.0"
  z="2.54*2.5"
  deltaphi="2*3.1415926535897"
  aunit="rad"
  lunit="cm"/>
</solids>
<structure>
 <volume name="vol_PMT_TPBCoating">
  <materialref ref="TPB"/>
  <solidref ref="PMT_TPBCoating"/>
 </volume>
 <volume name="vol_PMT_AcrylicPlate">
  <materialref ref="Acrylic"/>
  <solidref ref="PMT_AcrylicPlate"/>
 </volume>
 <volume name="vol_PMT_Stalk">
  <materialref ref="Glass"/>
  <solidref ref="PMT_Stalk"/>
 </volume>
 <volume name="vol_PMT_SteelBase">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="PMT_SteelBase"/>
 </volume>
 <volume name="vol_PMT_Underside">
  <materialref ref="Glass"/>
  <solidref ref="PMT_Underside"/>
 </volume>
 <volume name="vol_PMTSensitive">
  <materialref ref="LAr"/>
  <solidref ref="PMT_Lens"/>
 </volume>
 <volume name="volPMT">
  <materialref ref="LAr"/>
  <solidref ref="PMTVolume"/>
  <physvol>
   <volumeref ref="vol_PMT_TPBCoating"/>
   <position name="pos_PMT_TPBCoating" unit="cm" x="0" y="0" z="(5.5 * 2.54) - (0.5 * 0.005)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_AcrylicPlate"/>
   <position name="pos_PMT_AcrylicPlate" unit="cm" x="0" y="0" z="(5.5 * 2.54) - 0.01 - (0.5 * 0.2)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Stalk"/>
   <position name="pos_PMT_Stalk" unit="cm" x="0" y="0" z="(3.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_SteelBase"/>
   <position name="pos_PMT_SteelBase" unit="cm" x="0" y="0" z="(0.75 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMTSensitive"/>
   <position name="pos_PMT_Lens" unit="cm" x="0" y="0" z="(7.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Underside"/>
   <position name="pos_PMT_Underside" unit="cm" x="0" y="0" z="(7.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
 </volume>
</structure>
EOF

}




#Parameterize the steel cryostat that encloses the TPC.
sub gen_cryostat()
{
    # Set up the output file.
    $CRYOSTAT = "microboone/micro-cryostat" . $suffix . ".gdml";
    push (@gdmlFiles, $CRYOSTAT); # Add file to list of GDML fragments
    $CRYOSTAT = ">" . $CRYOSTAT;
    open(CRYOSTAT) or die("Could not open file $CRYOSTAT for writing");

    print CRYOSTAT <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <tube name="Cryostat" 
  rmax="$CryostatOuterRadius"
  z="$CryostatLength+2*$CryostatEndcapThickness"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
<tube name="SteelTube"
  rmin="$CryostatInnerRadius"
  rmax="$CryostatOuterRadius"
  z="$CryostatLength"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
EOF


	print CRYOSTAT <<EOF;
</solids>

<structure>
 <volume name="volSteelTube">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="SteelTube"/>
 </volume>
 <volume name="volCryostat">
  <materialref ref="LAr"/>
  <solidref ref="Cryostat"/>
  <physvol>
   <volumeref ref="volSteelTube"/>
   <position name="posSteelTube" unit="cm" x="0" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPC"/>
   <position name="posTPC" unit="cm" x="0.0" y="0.97" z="0"/>
  </physvol>
EOF


  @pmt_pos = ( ' x="-147.8"  y="3.21654"  z="-463.358"',
               ' x="-147.76" y="-52.6635" z="-403.121"',
               ' x="-147.8"  y="59.0965"  z="-403.121"',
               ' x="-147.76" y="-52.6635" z="-361.449"',
               ' x="-147.8"  y="59.0965"  z="-361.449"',
               ' x="-147.8"  y="3.21654"  z="-308.81"',
               ' x="-147.8"  y="3.21654"  z="-249.591"',
               ' x="-147.76" y="-52.6635" z="-194.487"',
               ' x="-147.8"  y="59.0965"  z="-194.487"',
               ' x="-147.76" y="-52.6635" z="-155.932"',
               ' x="-147.8"  y="59.0965"  z="-155.932"',
               ' x="-147.8"  y="3.21654"  z="-104.833"',
               ' x="-147.8"  y="3.21654"  z="-38.3029"',
               ' x="-147.76" y="-52.6635" z="13.6988"',
               ' x="-147.8"  y="59.0965"  z="13.6988"',
               ' x="-147.76" y="-52.6635" z="54.0326"',
               ' x="-147.8"  y="59.0965"  z="54.0326"',
               ' x="-147.8"  y="3.21654"  z="108.648"',
               ' x="-147.8"  y="3.21654"  z="175.178"',
               ' x="-147.76" y="-52.6635" z="229.136"',
               ' x="-147.8"  y="59.0965"  z="229.136"',
               ' x="-147.76" y="-52.6635" z="267.023"',
               ' x="-147.8"  y="59.0965"  z="267.023"',
               ' x="-147.8"  y="3.21654"  z="314.818"',
               ' x="-147.8"  y="3.21654"  z="385.004"',
               ' x="-147.76" y="-52.6635" z="434.169"',
               ' x="-147.8"  y="59.0965"  z="434.169"',
               ' x="-147.76" y="-52.6635" z="474.285"',
               ' x="-147.8"  y="59.0965"  z="474.285"',
               ' x="-147.8"  y="3.21654"  z="514.408"' );

  if ( $pmt_switch eq "on" ) {
    for ( $i=0; $i<30; ++$i ){
      print CRYOSTAT <<EOF;
  <physvol>
   <volumeref ref="volPMT"/>
   <position name="posPMT$i" unit="cm" @pmt_pos[$i]/>
   <rotationref ref="rPMTRotation1"/>
  </physvol>
EOF
    }
  }

	print CRYOSTAT <<EOF;
 </volume>
</structure>
</gdml>
EOF

   close(CRYOSTAT);
}


# Parameterize the cryostat's surroundings.
sub gen_enclosure()
{
    # Set up the output file.
    $GDML = "micro-enclosure" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="DetEnclosure" lunit="cm"
   x="kDetEnclosureWidth" y="$CryostatRadius" z="kDetEnclosureLength"
 />
</solids>

<structure>
 <volume name="volDetEnclosure">
  <materialref ref="Air"/>
  <solidref ref="DetEnclosure"/>
  <physvol>
   <volumeref ref="volCryostat"/>
   <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
  </physvol>
 </volume>
</structure>
</gdml>
EOF

   close(GDML);
}


sub gen_enclosure()
{
    # Set up the output file.
    $GDML = "micro-enclosure" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="DetEnclosure" lunit="cm"
   x="kDetEnclosureWidth" y="2*$CryostatOuterRadius+100" z="kDetEnclosureLength"
 />
</solids>

<structure>
 <volume name="volDetEnclosure">
  <materialref ref="Air"/>
  <solidref ref="DetEnclosure"/>
  <physvol>
   <volumeref ref="volCryostat"/>
   <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
  </physvol>
 </volume>
</structure>
</gdml>
EOF

   close(GDML);
}


# Parameterize the dirt mound that surrounds the enclosure.
sub gen_world()
{
    # Set up the output file.
    $GDML = "micro-world" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
  <box name="World" 
    lunit="cm" 
    x="kWorldW" 
    y="kWorldH" 
    z="kWorldL"/>
  <tube name="Ground"
    rmin="620*2.54"
    rmax="((50*12)+620)*2.54"
    z="41*12*2.54"
    deltaphi="360" 
    lunit="cm"
    aunit="deg"/>
  <tube name="ConcreteEnclosure"
    rmin="584*2.54"
    rmax="620*2.54"
    z="38*12*2.54"
    deltaphi="360" 
    lunit="cm"
    aunit="deg"/>
  <tube name="ConcreteEnclosureBottom"
    rmin="0"
    rmax="620*2.54"
    z="36*2.54"
    deltaphi="360" 
    lunit="cm"
    aunit="deg"/>
  <tube name="Overburden"
    rmin="0"
    rmax="584*2.54"
    z="10*12*2.54"
    deltaphi="360" 
    lunit="cm"
    aunit="deg"/>
</solids>

<structure>
  <volume name="volGround" >
    <materialref ref="Dirt" />
    <solidref ref="Ground" />
  </volume>
  <volume name="volOverburden" >
    <materialref ref="Dirt" />
    <solidref ref="Overburden" />
  </volume>
  <volume name="volConcreteEnclosure" >
    <materialref ref="Concrete" />
    <solidref ref="ConcreteEnclosure" />
  </volume>
  <volume name="volConcreteEnclosureBottom" >
    <materialref ref="Concrete" />
    <solidref ref="ConcreteEnclosureBottom" />
  </volume>
  <volume name="volWorld" >
    <materialref ref="Air"/> 
    <solidref ref="World"/>
    <physvol>
      <volumeref ref="volConcreteEnclosure"/>
      <position name="posConcreteEnclosure" unit="cm" x="0.5*kTPCDepth" y="36*2.54/2" z="0.5*kTPCLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
      <volumeref ref="volConcreteEnclosureBottom"/>
      <position name="posConcreteEnclosureBottom" unit="cm" x="0.5*kTPCDepth" y="-38*12*2.54/2" z="0.5*kTPCLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
      <volumeref ref="volGround"/>
      <position name="posGround" unit="cm" x="0.5*kTPCDepth" y="0" z="0.5*kTPCLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
      <volumeref ref="volOverburden"/>
      <position name="posOverburden" unit="cm" x="0.5*kTPCDepth" y="(41-10)*12*2.54/2" z="0.5*kTPCLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
      <volumeref ref="volDetEnclosure"/>
      <position name="posDetEnclosure" unit="cm" x="0.5*kTPCDepth" y="0" z="0.5*kTPCLength"/>
    </physvol>
  </volume> 
</structure>
</gdml>
EOF

   close(GDML);
}



sub write_fragments()
{
    # The output file is a list of the GDML sub-files created by this
    # script.

    if ( ! defined $output )
    {
	$output = "-"; # write to STDOUT 
    }

    # Set up the output file.
    $OUTPUT = ">" . $output;
    open(OUTPUT) or die("Could not open file $OUTPUT");

    print OUTPUT <<EOF;
<?xml version='1.0'?>

<!-- Input to Geometry/gdml/make_gdml.pl; define the GDML fragments
     that will be zipped together to create a detector description. 
     -->

<config>

   <constantfiles>

      <!-- These files contain GDML <constant></constant>
           blocks. They are read in separately, so they can be
           interpreted into the remaining GDML. See make_gdml.pl for
           more information. 
	   -->
	   
EOF

    foreach $filename (@defFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </constantfiles>

   <gdmlfiles>

      <!-- The GDML file fragments to be zipped together. -->

EOF

    foreach $filename (@gdmlFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </gdmlfiles>

</config>
EOF

    close(OUTPUT);
}
