unit UnitEQGD;
// ������ ����� Basic Parameters � ANSYS Icepak 12.0.
// ������ �������� ������ ��������
// 6 ������� 2016 ����.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, ExtCtrls, StdCtrls, Vcl.AppEvnts;

type
  TEGDForm = class(TForm)
    PanelGlobal: TPanel;
    GBCFlZ: TGroupBox;
    CBFlow: TCheckBox;
    RadioGroupFlowRegime: TRadioGroup;
    GroupBoxTurbulentModel: TGroupBox;
    ComboBoxturbulentmodel: TComboBox;
    ButtonidFlow: TButton;
    BEditTurb: TButton;
    GBGravity: TGroupBox;
    Lgx: TLabel;
    Lgy: TLabel;
    Lgz: TLabel;
    lblgx: TLabel;
    lblgy: TLabel;
    lbldz: TLabel;
    Egx: TEdit;
    Egy: TEdit;
    Egz: TEdit;
    ApplicationEvents1: TApplicationEvents;
    CheckBoxStaticStructural: TCheckBox;
    GroupBoxTemperature: TGroupBox;
    ComboBoxTemperature: TComboBox;
    //procedure ButtonTempClick(Sender: TObject);
    procedure RadioGroupFlowRegimeClick(Sender: TObject);
    procedure CBFlowClick(Sender: TObject);
    procedure ButtonidFlowClick(Sender: TObject);
    procedure CBIdCurFLzoneClick(Sender: TObject);
    procedure ComboBoxturbulentmodelChange(Sender: TObject);
    procedure BEditTurbClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure ComboBoxTemperatureChange(Sender: TObject);



  private
    { Private declarations }
    procedure patchstring(var s : String);
  public
    { Public declarations }
  end;

var
  EGDForm: TEGDForm;

implementation
uses
     VisualUnit, UnitSmagorinsky;

{$R *.dfm}

// ���������� ���������� � ��� ����� �� ������ ��������� ����������������.
{*
procedure TEGDForm.ButtonTempClick(Sender: TObject);
begin
   Laplas.egddata.itemper:=ComboBoxTemperature.ItemIndex;
end;
*}

{*
// ������� �� ����� ����������� ����� FLUID ���.
// ������ ��� ������� �.�. 5 ������� 2016 ��� ����������������� ����������
// ����� � ���� � �������� ��� ���� ������.
procedure TEGDForm.BmaxFluidDomainClick(Sender: TObject);
var
   i : Integer;
begin
   // ��������� ���������� ������������ FLUID ���
   Laplas.egddata.imaxflD:=CBMaxFluidDomain.ItemIndex;
   // ��������� ��� ����������� ����������� ������ ��-�� ��������� ��������
   // ������������� �������.
   SetLength(Laplas.egddata.myflmod,CBMaxFluidDomain.ItemIndex);

   // ������������� ����
   for i:=0 to Laplas.egddata.imaxflD-1 do
   begin
      Laplas.egddata.myflmod[i].xc:=0.0;
      Laplas.egddata.myflmod[i].yc:=0.0;
      Laplas.egddata.myflmod[i].zc:=0.0;
      Laplas.egddata.myflmod[i].iflow:=0; // �� ������� �������
      Laplas.egddata.myflmod[i].iflowregime:=0; // laminar
      Laplas.egddata.myflmod[i].iturbmodel:=0; // 0-ZEM, 1-Smagorinsky, 2 - RNG LES.
      // ������ �������������
      BEditTurb.Visible:=false; // ������ ������ ������������� ������� ������ ������������� ���������, �.�. ������� ZEM
      Laplas.egddata.myflmod[i].SmagConst:=0.151; // ��� Ck==1.8   (Ck ��������������� ��������� �����������).
      Laplas.egddata.myflmod[i].bSmagorinsky_Lilly:=true; // ������ �������������-�����.
      Laplas.egddata.myflmod[i].bfdelta:=true; // ���� ��������������� �����
      Laplas.egddata.myflmod[i].bsurface_roughness:=false; // �� ��������� ������������� ������.
      Laplas.egddata.myflmod[i].ipowerroughness:=2; // ���������� ������� � ������.
      Laplas.egddata.myflmod[i].roughness:=10.0; // micron (������������� ������ � ���).
      Laplas.egddata.myflmod[i].bSwirlamendment:=false;
      Laplas.egddata.myflmod[i].bSelectiveSmagorinsky:=False; // ������������� ������ �������������.
      //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].itypefiltr:=2; // ������ ��������.
      //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rSelectiveAngle:=15.0;
       Laplas.egddata.myflmod[0].itypefiltr:=2; // ������ ��������.
       Laplas.egddata.myflmod[0].rSelectiveAngle:=15.0;
   end;
   Application.MessageBox('information on the hydrodynamic DOMAIN cleared','Attantion!');

   if (CBMaxFluidDomain.ItemIndex=0) then
   begin
      GBCFlZ.Visible:=false;
   end
    else
   begin
      GBCFlZ.Visible:=true;

     // CBIdCurFLzone.Clear;
      //for i:=0 to Laplas.egddata.imaxflD-1 do
      //begin
        // CBIdCurFLzone.Items.Append(IntToStr(i+1));
      //end;
      //CBIdCurFLzone.ItemIndex:=0; // ������������� ������ FLUID ����

      // ���������� ��������� ������� �����
      //EditXC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].xc);
      //EditYC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].yc);
      //EditZC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].zc);

      // ����������� ��� �� ����������� ��������� �������������
     // if (Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow=1) then
      if (Laplas.egddata.myflmod[0].iflow=1) then
      begin
         CBFlow.Checked:=true;
      end
       else
      begin
         CBFlow.Checked:=false;
      end;

      // ����� �������
      //RadioGroupFlowRegime.ItemIndex:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflowregime;
      // ������ ��������������
      //ComboBoxturbulentmodel.ItemIndex:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iturbmodel;
      // ����� �������
      RadioGroupFlowRegime.ItemIndex:=Laplas.egddata.myflmod[0].iflowregime;
      // ������ ��������������
      ComboBoxturbulentmodel.ItemIndex:=Laplas.egddata.myflmod[0].iturbmodel;

   end;
end;
*}


procedure TEGDForm.RadioGroupFlowRegimeClick(Sender: TObject);
begin
   // ����� ������ �������:
   // ���������� ��� ������������.
   if (RadioGroupFlowRegime.ItemIndex=0) then
   begin
      // LAMINAR
      GroupBoxTurbulentModel.Visible:=false;
   end
    else
   begin
      // Turbulent
      GroupBoxTurbulentModel.Visible:=true;
   end;
end;

procedure TEGDForm.CBFlowClick(Sender: TObject);
begin
   // ������ ��� �� ������ ��������� �������������
   if (CBFlow.Checked) then
   begin
      // ��������� ������������� ��������.
      if ((ComboBoxTemperature.ItemIndex=2) or
       (ComboBoxTemperature.ItemIndex=3)) then
      begin
         ComboBoxTemperature.ItemIndex:=1;
      end;

      RadioGroupFlowRegime.Visible:=true;
      GBGravity.Visible:=true;
      // ���������� ��� ������������.
      if (RadioGroupFlowRegime.ItemIndex=0) then
      begin
         // LAMINAR
         GroupBoxTurbulentModel.Visible:=false;
      end
       else
      begin
         // Turbulent
         GroupBoxTurbulentModel.Visible:=true;
      end;
   end
    else
   begin
      // ��������� ������������� �� ��������.
      RadioGroupFlowRegime.Visible:=false;
      GroupBoxTurbulentModel.Visible:=false;
      GBGravity.Visible:=false;
      if ((ComboBoxTemperature.ItemIndex=0) and
       (CheckBoxStaticStructural.Checked=false)) then
      begin
         // �� ������������ �������������
         // �� ������������ ��������.
         // ������������ ����������� �������
         // ������������ ������.
         ComboBoxTemperature.ItemIndex:=1;
      end;

   end;
   // ���������� � ��� ������������ ��������� ������������� ��� ���.
   {*
   if (CBFlow.Checked) then
   begin
      // ��������� ��������� �������������.
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow:=1;
   end
    else
   begin
      // �� ��������� ��������� �������������.
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow:=0;
   end;
   *}
    if (CBFlow.Checked) then
   begin
      // ��������� ��������� �������������.
      Laplas.egddata.myflmod[0].iflow:=1;
      // ��������� ����� ����� �������� �����������������
      // ��������� ��������.
      Height:=488;
      ButtonidFlow.Top:=421;
      GBCFlZ.Height:=449;
      PanelGlobal.Height:=449;
   end
    else
   begin
      // �� ��������� ��������� �������������.
      Laplas.egddata.myflmod[0].iflow:=0;
      // �������� �������� � ������ ������ ����������
      // � ����������������� ��������������� �
      // ������� ����� ��� ������.
      GBCFlZ.Height:=135;
      Height:=136;
      PanelGlobal.Height:=136;
      ButtonidFlow.Top:=105;
   end;
end;

procedure TEGDForm.patchstring(var s : String);
var
   i : Integer;
begin
if (FormatSettings.DecimalSeparator='.') then
begin
   for i:=1 to length(s) do
   begin
      if (s[i]=',') then s[i]:='.';
   end;
end;

if (FormatSettings.DecimalSeparator=',') then
begin
   for i:=1 to length(s) do
   begin
      if (s[i]='.') then s[i]:=',';
   end;
end;
end;

procedure TEGDForm.ButtonidFlowClick(Sender: TObject);
var
  s1 : String;
  k : Integer;
   s : String;
  bOk : Boolean;
  c : Real;
  code : Integer;
  sforval : String;
begin
   // ������ ���������� ������� ���� �������.

   // ����������� ����������� �����������.
   s:=Trim(Egx.Text);
   patchstring(s);
   Egx.Text:=s;
   s:=Trim(Egy.Text);
   patchstring(s);
   Egy.Text:=s;
   s:=Trim(Egz.Text);
   patchstring(s);
   Egz.Text:=s;


   bOk:=False;
   //val(Trim(Egx.Text),c,code);
   sforval:='';
  sforval:=StringReplace(Egx.Text,',','.',[rfReplaceAll]);
   val(sforval,c,code);
   if (code=0) then
   begin
       //val(Trim(Egy.Text),c,code);
       sforval:='';
       sforval:=StringReplace(Egy.Text,',','.',[rfReplaceAll]);
       val(sforval,c,code);
      if (code=0) then
      begin
         //val(Trim(Egz.Text),c,code);
         sforval:='';
         sforval:=StringReplace(Egz.Text,',','.',[rfReplaceAll]);
         val(sforval,c,code);
         if (code=0) then
         begin
            bOk:=True;
         end;
      end;
   end;
   if (bOk) then
   begin
      Laplas.gx:=StrToFloat(Trim(Egx.Text));
      Laplas.gy:=StrToFloat(Trim(Egy.Text));
      Laplas.gz:=StrToFloat(Trim(Egz.Text));
   end
    else
   begin
      if (FormatSettings.DecimalSeparator='.') then
      begin
         Egx.Text:='0.0';
         Egy.Text:='0.0';
         Egz.Text:='0.0';
      end;
      if (FormatSettings.DecimalSeparator=',') then
      begin
         Egx.Text:='0,0';
         Egy.Text:='0,0';
         Egz.Text:='0,0';
      end;
      ShowMessage('Input Error, please repeat...');
   end;

    // ����� ���������� ������� ���� �������.


   // ���������� ���������� � ������� FLUID ����
   {*
   s1:=Trim(EditXC.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   EditXC.Text:=s1;

   s1:=Trim(EditYC.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   EditYC.Text:=s1;

   s1:=Trim(EditZC.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   EditZC.Text:=s1;


   // ���������� ������� �����
   Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].xc:=StrToFloat(EditXC.Text);
   Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].yc:=StrToFloat(EditYC.Text);
   Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].zc:=StrToFloat(EditZC.Text);
    *}

    // ��������� ��� ������� ��������� �������������.
    // ���������� ���������� � ��� ����� �� ������ ��������� ����������������.
    Laplas.egddata.itemper:=ComboBoxTemperature.ItemIndex;

   // ��������� ��� ������� ��������� ������ ���������.
   // ���������� ���������� � ��� ����� �� ������ �������� ��������� ���������.
   if (CheckBoxStaticStructural.Checked) then
   begin
      Laplas.egddata.iStaticStructural:=1;
   end
    else
   begin
      Laplas.egddata.iStaticStructural:=0;
   end;


    // � 5 ������� 2016 ���� � ��� ����� ���� ����������
    // FLUID ����, ��� ��������� ������� ������������
    // ����������������� ����������.
    // ���������� ������� �����
   Laplas.egddata.myflmod[0].xc:=0.0;
   Laplas.egddata.myflmod[0].yc:=0.0;
   Laplas.egddata.myflmod[0].zc:=0.0;

   {*
   // ���������� � ��� ������������ ��������� ������������� ��� ���.
   if (CBFlow.Checked) then
   begin
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow:=1;
   end
    else
   begin
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow:=0;
   end;

   // ����� ������� 0-���������� ��� 1-������������
   if (RadioGroupFlowRegime.ItemIndex=0) then
   begin
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflowregime:=0;
   end
    else
   begin
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflowregime:=1;
   end;

   // ���������� ���������� � ������ ��������������.
   Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iturbmodel:=ComboBoxturbulentmodel.ItemIndex;
   *}
   // ���������� � ��� ������������ ��������� ������������� ��� ���.
   if (CBFlow.Checked) then
   begin
      Laplas.egddata.myflmod[0].iflow:=1;
   end
    else
   begin
      Laplas.egddata.myflmod[0].iflow:=0;
   end;

   // ����� ������� 0-���������� ��� 1-������������
   if (RadioGroupFlowRegime.ItemIndex=0) then
   begin
      Laplas.egddata.myflmod[0].iflowregime:=0;
   end
    else
   begin
      Laplas.egddata.myflmod[0].iflowregime:=1;
   end;

   // ���������� ���������� � ������ ��������������.
   Laplas.egddata.myflmod[0].iturbmodel:=ComboBoxturbulentmodel.ItemIndex;

   Close();

end;

procedure TEGDForm.CBIdCurFLzoneClick(Sender: TObject);
begin
   // ����� ������� FLUID ����

   // ���������� ��������� ������� �����
  // EditXC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].xc);
   //EditYC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].yc);
   //EditZC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].zc);

   // ����������� ��� �� ����������� ��������� �������������
   //if (Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow=1) then
   if (Laplas.egddata.myflmod[0].iflow=1) then
   begin
      CBFlow.Checked:=true;
   end
    else
   begin
      CBFlow.Checked:=false;
   end;

   // ����� �������
   //RadioGroupFlowRegime.ItemIndex:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflowregime;
   // ������ ��������������
   //ComboBoxturbulentmodel.ItemIndex:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iturbmodel;
    // ����� �������
   RadioGroupFlowRegime.ItemIndex:=Laplas.egddata.myflmod[0].iflowregime;
   // ������ ��������������
   ComboBoxturbulentmodel.ItemIndex:=Laplas.egddata.myflmod[0].iturbmodel;
end;

// ������ ������� ��� ��������� ������ ��������������
// ������� ������ �������������.
procedure TEGDForm.ComboBoxTemperatureChange(Sender: TObject);
begin
   if ((ComboBoxTemperature.ItemIndex=0)
   and(CBFlow.Checked=false)and
   (CheckBoxStaticStructural.Checked=false)) then
   begin
       // ������ �� ��������� ������ ������� ������������ ������.
       ComboBoxTemperature.ItemIndex:=1; // �� ��� ������� �����������,
       // ���� ��������� ������������� � ��������.
   end;
   if ((ComboBoxTemperature.ItemIndex=2)and
   (CBFlow.Checked=true)and
   (CheckBoxStaticStructural.Checked=false)) then
   begin
      // Finite Element Method -> Finite Volume Method.
      ComboBoxTemperature.ItemIndex:=1;
   end;
    if ((ComboBoxTemperature.ItemIndex=3)and
   (CBFlow.Checked=true)and
   (CheckBoxStaticStructural.Checked=false)) then
   begin
      // Network T Solver -> Finite Volume Method.
      ComboBoxTemperature.ItemIndex:=1;
   end;
end;

procedure TEGDForm.ComboBoxturbulentmodelChange(Sender: TObject);
begin
   // ����� ������ ��������������
   case ComboBoxturbulentmodel.ItemIndex of
     0 : begin
            // ZEM (RANS)
            BEditTurb.Visible:=false;
         end;
     1 : begin
            // Smagorinsky
            BEditTurb.Visible:=true;
         end;
     2 : begin
            // RNG LES
            BEditTurb.Visible:=false;
         end;
     3 : begin
            // Spalart Allmares (RANS) [1992]
            BEditTurb.Visible:=false;
         end;
     4 : begin
            // SST ������� (RANS) [1993]
            BEditTurb.Visible:=false;
         end;
     5 : begin
            // ����������� ������ �� ������
            // Standart K-Epsilon (RANS) [2001]
            BEditTurb.Visible:=false;
         end;
     6: begin
           // ������ ������� � ��������
           // ��������� ������������� ��������  [2009]
           BEditTurb.Visible:=false;
         end;
   end;
end;




// ������������� ����� GravityForm ��� �������� �����.
procedure TEGDForm.FormCreate(Sender: TObject);
begin
    Egx.Text:=FloatToStr(Laplas.gx);
    Egy.Text:=FloatToStr(Laplas.gy);
    Egz.Text:=FloatToStr(Laplas.gz);
end;

// ������ ����� �������������.
procedure TEGDForm.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

procedure TEGDForm.BEditTurbClick(Sender: TObject);
begin
{*
   // ��������� ������� ��� ��� 5 ������� 2016 �� ����� ��� ����������������� ����������
   // � ���� ������������  � ��������� ��� ���� ��������.
   // ������� �� ������ �������������� ������ �������������.
   if (ComboBoxturbulentmodel.ItemIndex=1) then
   begin
      if (Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iDynamicStressGermano=1) then
      begin
         // ������ ������� ��������.
         FormSmagorinsky.PanelCs.Visible:=false;
         FormSmagorinsky.PanelLimitersCs.Visible:=true;
         FormSmagorinsky.CBDynamicStress.Checked:=true;
         if (Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].iLimitersCs=0) then
         begin
            FormSmagorinsky.CBLimitersCS.Checked:=false;
            FormSmagorinsky.Panel_user_limiters.Visible:=false;
         end
          else
         begin
            FormSmagorinsky.CBLimitersCS.Checked:=true;
            // ������� Cs ��������
            FormSmagorinsky.Panel_user_limiters.Visible:=true;
            FormSmagorinsky.EminCs.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].minCs);
            FormSmagorinsky.EmaxCs.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].maxCs);
         end;
      end
       else
      begin
         // ������ ������� ���������.
         FormSmagorinsky.PanelCs.Visible:=true;
         FormSmagorinsky.Esmagconst.Text:=FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].SmagConst);
         FormSmagorinsky.PanelLimitersCs.Visible:=false;
         FormSmagorinsky.CBDynamicStress.Checked:=false;
      end;

      FormSmagorinsky.CBcorrectnonuniformgrid.Checked:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].bfdelta;
      FormSmagorinsky.CBSmagLilly.Checked:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].bSmagorinsky_Lilly;
      FormSmagorinsky.CBsurfaceroughness.Checked:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].bsurface_roughness;
      if (FormSmagorinsky.CBsurfaceroughness.Checked) then
      begin
         // ����� ������������ ����� ��� ������� ������������� ������.
         FormSmagorinsky.Panrougness.Visible:=true;
         FormSmagorinsky.Eroughnessmicron.Text:=FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].roughness);
         case Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].ipowerroughness of
           1 : begin
                  FormSmagorinsky.CBpowermodel.ItemIndex:=0;
               end;
           2 : begin
                  FormSmagorinsky.CBpowermodel.ItemIndex:=1;
               end;
         end;
      end
      else
      begin
         FormSmagorinsky.Panrougness.Visible:=false;
      end;
      // �������� ��� ������� � ��������� ����� ����.
      FormSmagorinsky.CBSwirlAmendment.Checked:=Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bSwirlamendment;
      if (FormSmagorinsky.CBSwirlAmendment.Checked) then
      begin
         FormSmagorinsky.PanRichardson.Visible:=true;
         FormSmagorinsky.ERimult.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rRimult);
      end
        else
      begin
         FormSmagorinsky.PanRichardson.Visible:=false;
      end;
      // ��������� ������ Selective Smagorinsky
      FormSmagorinsky.CBSelectiveSmagorinsky.Checked:=Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bSelectiveSmagorinsky;
      if (FormSmagorinsky.CBSelectiveSmagorinsky.Checked) then
      begin
         FormSmagorinsky.PSelectiveSmag.Visible:=true;
         FormSmagorinsky.RGtypefiltr.ItemIndex:=Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].itypefiltr;
         FormSmagorinsky.Eangle.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rSelectiveAngle);
      end
       else
      begin
         FormSmagorinsky.PSelectiveSmag.Visible:=false;
      end;
      FormSmagorinsky.ShowModal;
   end;
   *}
   // ������� �� ������ �������������� ������ �������������.
   if (ComboBoxturbulentmodel.ItemIndex=1) then
   begin
      if (Laplas.egddata.myflmod[0].iDynamicStressGermano=1) then
      begin
         // ������ ������� ��������.
         FormSmagorinsky.PanelCs.Visible:=false;
         FormSmagorinsky.PanelLimitersCs.Visible:=true;
         FormSmagorinsky.CBDynamicStress.Checked:=true;
         if (Laplas.egddata.myflmod[0].iLimitersCs=0) then
         begin
            FormSmagorinsky.CBLimitersCS.Checked:=false;
            FormSmagorinsky.Panel_user_limiters.Visible:=false;
         end
          else
         begin
            FormSmagorinsky.CBLimitersCS.Checked:=true;
            // ������� Cs ��������
            FormSmagorinsky.Panel_user_limiters.Visible:=true;
            FormSmagorinsky.EminCs.Text:=FloatToStr(Laplas.egddata.myflmod[0].minCs);
            FormSmagorinsky.EmaxCs.Text:=FloatToStr(Laplas.egddata.myflmod[0].maxCs);
         end;
      end
       else
      begin
         // ������ ������� ���������.
         FormSmagorinsky.PanelCs.Visible:=true;
         FormSmagorinsky.Esmagconst.Text:=FloatToStr(Laplas.egddata.myflmod[0].SmagConst);
         FormSmagorinsky.PanelLimitersCs.Visible:=false;
         FormSmagorinsky.CBDynamicStress.Checked:=false;
      end;

      FormSmagorinsky.CBcorrectnonuniformgrid.Checked:=Laplas.egddata.myflmod[0].bfdelta;
      FormSmagorinsky.CBSmagLilly.Checked:=Laplas.egddata.myflmod[0].bSmagorinsky_Lilly;
      FormSmagorinsky.CBsurfaceroughness.Checked:=Laplas.egddata.myflmod[0].bsurface_roughness;
      if (FormSmagorinsky.CBsurfaceroughness.Checked) then
      begin
         // ����� ������������ ����� ��� ������� ������������� ������.
         FormSmagorinsky.Panrougness.Visible:=true;
         FormSmagorinsky.Eroughnessmicron.Text:=FloatToStr(Laplas.egddata.myflmod[0].roughness);
         case Laplas.egddata.myflmod[0].ipowerroughness of
           1 : begin
                  FormSmagorinsky.CBpowermodel.ItemIndex:=0;
               end;
           2 : begin
                  FormSmagorinsky.CBpowermodel.ItemIndex:=1;
               end;
         end;
      end
      else
      begin
         FormSmagorinsky.Panrougness.Visible:=false;
      end;
      // �������� ��� ������� � ��������� ����� ����.
      FormSmagorinsky.CBSwirlAmendment.Checked:=Laplas.egddata.myflmod[0].bSwirlamendment;
      if (FormSmagorinsky.CBSwirlAmendment.Checked) then
      begin
         FormSmagorinsky.PanRichardson.Visible:=true;
         FormSmagorinsky.ERimult.Text:=FloatToStr(Laplas.egddata.myflmod[0].rRimult);
      end
        else
      begin
         FormSmagorinsky.PanRichardson.Visible:=false;
      end;
      // ��������� ������ Selective Smagorinsky
      FormSmagorinsky.CBSelectiveSmagorinsky.Checked:=Laplas.egddata.myflmod[0].bSelectiveSmagorinsky;
      if (FormSmagorinsky.CBSelectiveSmagorinsky.Checked) then
      begin
         FormSmagorinsky.PSelectiveSmag.Visible:=true;
         FormSmagorinsky.RGtypefiltr.ItemIndex:=Laplas.egddata.myflmod[0].itypefiltr;
         FormSmagorinsky.Eangle.Text:=FloatToStr(Laplas.egddata.myflmod[0].rSelectiveAngle);
      end
       else
      begin
         FormSmagorinsky.PSelectiveSmag.Visible:=false;
      end;
      FormSmagorinsky.ShowModal;
   end;
end;

end.
