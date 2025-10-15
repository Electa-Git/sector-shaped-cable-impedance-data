// File created with GetDP.jl: https://github.com/Electa-Git/GetDP.jl.

DefineConstant[
active_con = {1, Choices{1,9999}, Name "Input/Active conductor", Visible 1}];

Group{
  DomainInf = Region[ {300100201,  300200102} ];// Domain transformation to infinity
  Con_1 = Region[ {100101103} ];// cable_1_core1_con_aluminum
  Con_2 = Region[ {100102103} ];// cable_1_core2_con_aluminum
  Con_3 = Region[ {100103103} ];// cable_1_core3_con_aluminum
  Con_4 = Region[ {100104104} ];// cable_1_neutral_con_copper
  Conductors = Region[ {100101103,  100102103,  100103103,  100104104} ];
  DomainC = Region[ {} ];// All conductor materials
  DomainCC = Region[ {} ];// All non-conductor materials
  DomainActive = Region[ {Con~{active_con}} ];// Sources
  DomainInactive = Region[ {Conductors - Con~{active_con}} ];// Conductors set to zero energization
  DomainC += Region[ {100101103} ];// cable_1_core1_con_aluminum
  DomainC += Region[ {100102103} ];// cable_1_core2_con_aluminum
  DomainC += Region[ {100103103} ];// cable_1_core3_con_aluminum
  DomainC += Region[ {100104104} ];// cable_1_neutral_con_copper
  DomainC += Region[ {200200102} ];// layer_2_earth_con_material_rho=100.0_epsr=1.0_mu=1.0
  DomainC += Region[ {300200102} ];// infshell_2_earth_con_material_rho=100.0_epsr=1.0_mu=1.0
  DomainCC += Region[ {100101201} ];// cable_1_core1_ins_air
  DomainCC += Region[ {100102201} ];// cable_1_core2_ins_air
  DomainCC += Region[ {100103201} ];// cable_1_core3_ins_air
  DomainCC += Region[ {100104201} ];// cable_1_neutral_ins_air
  DomainCC += Region[ {200100201} ];// layer_1_air_ins_air
  DomainCC += Region[ {300100201} ];// infshell_1_air_ins_air
  Domain_Ele = Region[ {DomainCC,  DomainC} ];
  Sur_Dirichlet_Ele = Region[ {1002201,  1001201} ];

}
Function{
  // Material properties for region 100104104: cable_1_neutral_con_copper
  nu[Region[{100104104}]] = 795774.7154594767;
  sigma[Region[{100104104}]] = 5.312815223651904e7;
  epsilon[Region[{100104104}]] = 8.8541878128e-12;
  // Material properties for region 100101103: cable_1_core1_con_aluminum
  nu[Region[{100101103}]] = 795774.7154594767;
  sigma[Region[{100101103}]] = 3.369311073797809e7;
  epsilon[Region[{100101103}]] = 8.8541878128e-12;
  /*  Material properties for region 200200102:
  layer_2_earth_con_material_rho=100.0_epsr=1.0_mu=1.0 */
  nu[Region[{200200102}]] = 795774.7154594767;
  sigma[Region[{200200102}]] = 0.01;
  epsilon[Region[{200200102}]] = 8.8541878128e-12;
  // Material properties for region 100102201: cable_1_core2_ins_air
  nu[Region[{100102201}]] = 795774.7154594767;
  sigma[Region[{100102201}]] = 0.0;
  epsilon[Region[{100102201}]] = 7.08335025024e-11;
  // Material properties for region 300100201: infshell_1_air_ins_air
  nu[Region[{300100201}]] = 795774.7154594767;
  sigma[Region[{300100201}]] = 0.0;
  epsilon[Region[{300100201}]] = 8.8541878128e-12;
  // Material properties for region 100101201: cable_1_core1_ins_air
  nu[Region[{100101201}]] = 795774.7154594767;
  sigma[Region[{100101201}]] = 0.0;
  epsilon[Region[{100101201}]] = 7.08335025024e-11;
  // Material properties for region 100104201: cable_1_neutral_ins_air
  nu[Region[{100104201}]] = 795774.7154594767;
  sigma[Region[{100104201}]] = 0.0;
  epsilon[Region[{100104201}]] = 8.8541878128e-12;
  // Material properties for region 100103103: cable_1_core3_con_aluminum
  nu[Region[{100103103}]] = 795774.7154594767;
  sigma[Region[{100103103}]] = 3.369311073797809e7;
  epsilon[Region[{100103103}]] = 8.8541878128e-12;
  // Material properties for region 100102103: cable_1_core2_con_aluminum
  nu[Region[{100102103}]] = 795774.7154594767;
  sigma[Region[{100102103}]] = 3.369311073797809e7;
  epsilon[Region[{100102103}]] = 8.8541878128e-12;
  // Material properties for region 100103201: cable_1_core3_ins_air
  nu[Region[{100103201}]] = 795774.7154594767;
  sigma[Region[{100103201}]] = 0.0;
  epsilon[Region[{100103201}]] = 7.08335025024e-11;
  // Material properties for region 200100201: layer_1_air_ins_air
  nu[Region[{200100201}]] = 795774.7154594767;
  sigma[Region[{200100201}]] = 0.0;
  epsilon[Region[{200100201}]] = 8.8541878128e-12;
  /*  Material properties for region 300200102:
  infshell_2_earth_con_material_rho=100.0_epsr=1.0_mu=1.0 */
  nu[Region[{300200102}]] = 795774.7154594767;
  sigma[Region[{300200102}]] = 0.01;
  epsilon[Region[{300200102}]] = 8.8541878128e-12;

}
Function{
  Freq = 50.0;
  UnitAmplitude = 1.0;

}
Constraint{
  { Name ScalarPotential_2D; Type Assign;
    Case {
      { Region DomainInactive; Value 0.0; }
      { Region Con~{active_con}; Value UnitAmplitude; }
      { Region Sur_Dirichlet_Ele; Value 0.0; }
    }
  }
  { Name Charge_2D; Type Assign;
    Case {
    }
  }
}
FunctionSpace {
  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node; Support Domain_Ele; Entity NodesOf[ All, Not Conductors ]; }
      { Name sf; NameOfCoef vf; Function BF_GroupOfNodes; Support Domain_Ele; Entity GroupsOfNodesOf[ Conductors ]; }
    }
    GlobalQuantity {
      { Name U; Type AliasOf; NameOfCoef vf;}
      { Name Q; Type AssociatedWith; NameOfCoef vf;}
    }
    Constraint {
      { NameOfCoef U; EntityType Region; NameOfConstraint ScalarPotential_2D;}
      { NameOfCoef Q; EntityType Region; NameOfConstraint Charge_2D;}
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint ScalarPotential_2D;}
    }
  }
}

Jacobian{
  { Name Vol; Case {
   {  Region DomainInf; Jacobian VolSphShell{503.2921210448703, 629.1151513060879, 0.0, 0.0, 0.0}; } 
   {  Region All; Jacobian Vol; } 
  } }
  { Name Sur; Case {
   {  Region All; Jacobian Sur; } 
  } }
}
Integration {
  { Name I1;
    Case {
      { Type Gauss;
      Case {
         {  GeoElement Point; NumberOfPoints 1; } 
         {  GeoElement Line; NumberOfPoints 4; } 
         {  GeoElement Triangle; NumberOfPoints 4; } 
         {  GeoElement Quadrangle; NumberOfPoints 4; } 
      }
      }
    }
  }
}
Formulation {
  { Name Electrodynamics_v; Type FemEquation;
      Quantity {
        { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
        { Name U; Type Global; NameOfSpace Hgrad_v_Ele [U]; }
        { Name Q; Type Global; NameOfSpace Hgrad_v_Ele [Q]; }
      }
      Equation {
        Galerkin { [ sigma[] * Dof{d v} , {d v} ];
            In Domain_Ele;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { DtDof[ epsilon[] * Dof{d v} , {d v} ];
            In DomainCC;
            Jacobian Vol;
            Integration I1;
        }
        GlobalTerm {  [ Dof{Q} , {U} ] ;
            In Conductors;
        }
      }
  }
}
Resolution {
  { Name Electrodynamics;
    System {
      { Name Sys_Ele; NameOfFormulation Electrodynamics_v; Type Complex; Frequency Freq; }
    }
    Operation {
      CreateDir["results/electrodynamics"];
      Generate[Sys_Ele];
      Solve[Sys_Ele];
      SaveSolution[Sys_Ele];
      PostOperation[LineParams];
    }
  }
}

PostProcessing {
  { Name EleDyn_v; NameOfFormulation Electrodynamics_v;
      PostQuantity {
      { Name v; Value {
          Term { [ {v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name e; Value {
          Term { [ -{d v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name em; Value {
          Term { [ Norm[-{d v}] ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name d; Value {
          Term { [ -epsilon[] * {d v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name dm; Value {
          Term { [ Norm[-epsilon[] * {d v}] ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name j; Value {
          Term { [ -sigma[] * {d v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name jm; Value {
          Term { [ Norm[-sigma[] * {d v}] ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name jtot; Value {
          Term { Type Global; [ -sigma[] * {d v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
          Term { Type Global; [ -epsilon[] * Dt[{d v}] ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name U; Value {
          Term { [ {U} ];
              In Domain_Ele;
          }
      }}
      { Name Q; Value {
          Term { [ {Q} ];
              In Domain_Ele;
          }
      }}
      { Name Y; Value {
          Term { [ -{Q} ];
              In Domain_Ele;
          }
      }}
      }
  }
}

PostOperation {
  { Name Field_Maps; NameOfPostProcessing EleDyn_v;
      Operation {
          Print[ v, OnElementsOf Domain_Ele, File StrCat[ "results/electrodynamics\v_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ em, OnElementsOf Domain_Ele, Name "|E| [V/m]", File StrCat[ "results/electrodynamics\em_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ dm, OnElementsOf Domain_Ele, Name "|D| [A/m²]", File StrCat[ "results/electrodynamics\dm_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ e, OnElementsOf Domain_Ele, Name "E [V/m]", File StrCat[ "results/electrodynamics\e_", Sprintf("%g",active_con), ".pos" ] ];
      }
  }
  { Name LineParams; NameOfPostProcessing EleDyn_v;
      Operation {
          Print[ Y, OnRegion Conductors, Format Table, File "results/electrodynamics\Y.dat", AppendToExistingFile (active_con > 1 ? 1 : 0) ];
      }
  }
}
