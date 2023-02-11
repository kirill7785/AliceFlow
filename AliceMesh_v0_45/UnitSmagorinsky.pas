unit UnitSmagorinsky;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TFormSmagorinsky = class(TForm)
    Panelmain: TPanel;
    CBcorrectnonuniformgrid: TCheckBox;
    CBSmagLilly: TCheckBox;
    GBsurfacerougnessparam: TGroupBox;
    BApply: TButton;
    GBSwirlDominatedFlow: TGroupBox;
    CBSwirlAmendment: TCheckBox;
    Panrougness: TPanel;
    LRoughness: TLabel;
    Eroughnessmicron: TEdit;
    Lroghness: TLabel;
    Lpower: TLabel;
    CBpowermodel: TComboBox;
    CBsurfaceroughness: TCheckBox;
    PanRichardson: TPanel;
    LRiMult: TLabel;
    ERimult: TEdit;
    GBSelectiveSmagorinsky: TGroupBox;
    CBSelectiveSmagorinsky: TCheckBox;
    PSelectiveSmag: TPanel;
    RGtypefiltr: TRadioGroup;
    LselectiveAngle: TLabel;
    Eangle: TEdit;
    Ldegree: TLabel;
    GBConstSmagorinsky: TGroupBox;
    PanelCs: TPanel;
    Lsmagconst: TLabel;
    Esmagconst: TEdit;
    CBDynamicStress: TCheckBox;
    PanelLimitersCs: TPanel;
    CBLimitersCs: TCheckBox;
    Panel_user_limiters: TPanel;
    LminCs: TLabel;
    LmaxCs: TLabel;
    RGtypefiltrGermanoModel: TRadioGroup;
    EminCs: TEdit;
    EmaxCs: TEdit;
    procedure BApplyClick(Sender: TObject);
    procedure CBsurfaceroughnessClick(Sender: TObject);
    procedure CBSwirlAmendmentClick(Sender: TObject);
    procedure CBSelectiveSmagorinskyClick(Sender: TObject);
    procedure CBDynamicStressClick(Sender: TObject);
    procedure CBLimitersCsClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormSmagorinsky: TFormSmagorinsky;

implementation
uses
     VisualUnit, UnitEQGD; // ���������� ������ �������� ������.
{$R *.dfm}

procedure TFormSmagorinsky.BApplyClick(Sender: TObject);
var
  s1 : String;
  k : Integer;
begin

    s1:=Trim(Esmagconst.Text);
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
   Esmagconst.Text:=s1;

    s1:=Trim(EminCs.Text);
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
   EminCs.Text:=s1;

    s1:=Trim(EmaxCs.Text);
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
   EmaxCs.Text:=s1;



    s1:=Trim(Eroughnessmicron.Text);
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
   Eroughnessmicron.Text:=s1;

    s1:=Trim(ERimult.Text);
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
   ERimult.Text:=s1;


    s1:=Trim(Eangle.Text);
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
   Eangle.Text:=s1;


    // ������� ������������� ������� ������.
    if (not(CBDynamicStress.Checked)) then
       begin
          if (length(Esmagconst.Text)=0) then
          begin
             //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].SmagConst:=0.1;
             Laplas.egddata.myflmod[0].SmagConst:=0.1;
             Laplas.MainMemo.Lines.Add('Smagorinsky const is equal 0.1');
             Application.MessageBox('Smagorinsky const is equal 0.1','Attantion!',MB_OK);
          end
           else
          begin
            // Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].SmagConst:=StrToFloat(Esmagconst.Text); // ���������� �������������.
             Laplas.egddata.myflmod[0].SmagConst:=StrToFloat(Esmagconst.Text); // ���������� �������������.
          end;
       end;
    if (CBDynamicStress.Checked) then
    begin
       // ������ ������� ��������.
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].iDynamicStressGermano:=1;
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].itypeFiltrGermano:=RGtypefiltrGermanoModel.ItemIndex; // ��� ������� ������� ������������ � ������ �������.
       Laplas.egddata.myflmod[0].iDynamicStressGermano:=1;
       Laplas.egddata.myflmod[0].itypeFiltrGermano:=RGtypefiltrGermanoModel.ItemIndex; // ��� ������� ������� ������������ � ������ �������.
       if (CBLimitersCs.Checked) then
       begin
          // ������������ ������������ �� ��������� �������������.
          //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].iLimitersCs:=1;
          // ���� �������������� ���������� ������������� ��������. ��������� ��������
          // �������������.
          //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].minCs:=StrToFloat(EminCs.Text);
          //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].maxCs:=StrToFloat(EmaxCs.Text);
           Laplas.egddata.myflmod[0].iLimitersCs:=1;
          // ���� �������������� ���������� ������������� ��������. ��������� ��������
          // �������������.
          Laplas.egddata.myflmod[0].minCs:=StrToFloat(EminCs.Text);
          Laplas.egddata.myflmod[0].maxCs:=StrToFloat(EmaxCs.Text);
       end
        else
       begin
          // �������� ����������� �� ��������� ������������� �� ������������,
          // ������������� �������� ��������� ����� �������� � ������������ �������.
          // Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].iLimitersCs:=0;
          //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].minCs:=-1e20; // ������ ��������� �� �������� ����������������
          //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].maxCs:=1e23;
           Laplas.egddata.myflmod[0].iLimitersCs:=0;
          Laplas.egddata.myflmod[0].minCs:=-1e20; // ������ ��������� �� �������� ����������������
          Laplas.egddata.myflmod[0].maxCs:=1e23;
       end;
    end
     else
    begin
       // ������ ������� ���������. �� ������������ ��������� ��������� �������������.
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].iDynamicStressGermano:=0;
       // �� �� ���������� ������� ����������� �� ��������� ������������� � ������ �����
       // ������������ ������ ������� �� ������������.
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].iLimitersCs:=0;
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].minCs:=-1.0e20; // ������ ��������� �� �������� ����������������
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].maxCs:=1.0e23;
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].itypeFiltrGermano:=2; // ������ ��������.
       Laplas.egddata.myflmod[0].iDynamicStressGermano:=0;
       // �� �� ���������� ������� ����������� �� ��������� ������������� � ������ �����
       // ������������ ������ ������� �� ������������.
       Laplas.egddata.myflmod[0].iLimitersCs:=0;
       Laplas.egddata.myflmod[0].minCs:=-1.0e20; // ������ ��������� �� �������� ����������������
       Laplas.egddata.myflmod[0].maxCs:=1.0e23;
       Laplas.egddata.myflmod[0].itypeFiltrGermano:=2; // ������ ��������.

    end;
    //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bfdelta:=CBcorrectnonuniformgrid.Checked; // ���� ��������������� �����.
    //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bSmagorinsky_Lilly:=CBSmagLilly.Checked; // ������ �������������-�����.
    //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bsurface_roughness:=CBsurfaceroughness.Checked; // ���� ������������� ������.
    Laplas.egddata.myflmod[0].bfdelta:=CBcorrectnonuniformgrid.Checked; // ���� ��������������� �����.
    Laplas.egddata.myflmod[0].bSmagorinsky_Lilly:=CBSmagLilly.Checked; // ������ �������������-�����.
    Laplas.egddata.myflmod[0].bsurface_roughness:=CBsurfaceroughness.Checked; // ���� ������������� ������.


    if (CBsurfaceroughness.Checked) then
    begin
       // ����� ��������� ������������� ������.
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].roughness:=StrToFloat(Eroughnessmicron.Text); // ������������� � ���.
        Laplas.egddata.myflmod[0].roughness:=StrToFloat(Eroughnessmicron.Text); // ������������� � ���.
       case CBpowermodel.ItemIndex of
         0 : begin
                //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].ipowerroughness:=1;
                  Laplas.egddata.myflmod[0].ipowerroughness:=1;
             end;
         1 : begin
                //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].ipowerroughness:=2;
                 Laplas.egddata.myflmod[0].ipowerroughness:=2;
             end;
       end;
    end;
    // ���� ������� � ��������� ����� ����.
    //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bSwirlamendment:=CBSwirlAmendment.Checked;
     Laplas.egddata.myflmod[0].bSwirlamendment:=CBSwirlAmendment.Checked;
    if (CBSwirlAmendment.Checked) then
    begin
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rRimult:=StrToFloat(ERimult.Text);
        Laplas.egddata.myflmod[0].rRimult:=StrToFloat(ERimult.Text);
    end;
    // ��������� ��� ���������� ������ Selective Smagorinsky.
    //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bSelectiveSmagorinsky:=CBSelectiveSmagorinsky.Checked;
    Laplas.egddata.myflmod[0].bSelectiveSmagorinsky:=CBSelectiveSmagorinsky.Checked;
    if (CBSelectiveSmagorinsky.Checked) then
    begin
       // ��� ������������� �������.
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].itypefiltr:=RGtypefiltr.ItemIndex;
       //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rSelectiveAngle:=StrToFloat(Eangle.Text);
       Laplas.egddata.myflmod[0].itypefiltr:=RGtypefiltr.ItemIndex;
       Laplas.egddata.myflmod[0].rSelectiveAngle:=StrToFloat(Eangle.Text);
    end;
end;

// ������ ������ ����� �������������.
procedure TFormSmagorinsky.CBsurfaceroughnessClick(Sender: TObject);
begin
   // �������� ��� ������ ������ ��� ����� �������������.
   if (CBsurfaceroughness.Checked) then
   begin
      Panrougness.Visible:=true;
      //Eroughnessmicron.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].roughness);
      Eroughnessmicron.Text:=FloatToStr(Laplas.egddata.myflmod[0].roughness);
      {*
      case Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].ipowerroughness of
        1 : begin
               CBpowermodel.ItemIndex:=0;
            end;
        2 : begin
               CBpowermodel.ItemIndex:=1;
            end;
       end;
       *}

        case Laplas.egddata.myflmod[0].ipowerroughness of
        1 : begin
               CBpowermodel.ItemIndex:=0;
            end;
        2 : begin
               CBpowermodel.ItemIndex:=1;
            end;
       end;
   end
    else
   begin
      Panrougness.Visible:=false;
   end;
end;

// ������ ��� �������� ������ �� ���� ������� � ��������� ����� ����.
procedure TFormSmagorinsky.CBSwirlAmendmentClick(Sender: TObject);
begin
   // ������ ��� ��������� ������ ��� ����� ������� � ��������� ����� ����.
   if (CBSwirlAmendment.Checked) then
   begin
      PanRichardson.Visible:=true;
      //ERimult.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rRimult);
      ERimult.Text:=FloatToStr(Laplas.egddata.myflmod[0].rRimult);
   end
     else
   begin
      PanRichardson.Visible:=false;
   end;
end;

procedure TFormSmagorinsky.CBSelectiveSmagorinskyClick(Sender: TObject);
begin
   // ������ ��� �������� ������ ���������� ��� �����������
   // ������ Selective Smagorinsky.
   if (CBSelectiveSmagorinsky.Checked) then
   begin
      PSelectiveSmag.Visible:=true;
      //RGtypefiltr.ItemIndex:=Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].itypefiltr;
      //Eangle.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rSelectiveAngle);
      RGtypefiltr.ItemIndex:=Laplas.egddata.myflmod[0].itypefiltr;
      Eangle.Text:=FloatToStr(Laplas.egddata.myflmod[0].rSelectiveAngle);
   end
    else
   begin
      PSelectiveSmag.Visible:=false;
   end;
end;

// ���������� ��� ������������ ������������ ������.
procedure TFormSmagorinsky.CBDynamicStressClick(Sender: TObject);
begin
   // ��������� ��� ���������� ������������ ������
   if (CBDynamicStress.Checked) then
   begin
      // ������ ������� ��������.
      PanelCs.Visible:=false;
      PanelLimitersCs.Visible:=true;
      // if (Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].iLimitersCs=0) then
      if (Laplas.egddata.myflmod[0].iLimitersCs=0) then
      begin
         CBLimitersCS.Checked:=false;
         Panel_user_limiters.Visible:=false;
      end
       else
      begin
         CBLimitersCS.Checked:=true;
         // ������� Cs ��������
         Panel_user_limiters.Visible:=true;
         //EminCs.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].minCs);
         //EmaxCs.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].maxCs);
          EminCs.Text:=FloatToStr(Laplas.egddata.myflmod[0].minCs);
         EmaxCs.Text:=FloatToStr(Laplas.egddata.myflmod[0].maxCs);
      end;
   end
    else
   begin
      // ������ ������� ���������.
      PanelCs.Visible:=true;
      //Esmagconst.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].SmagConst);
       Esmagconst.Text:=FloatToStr(Laplas.egddata.myflmod[0].SmagConst);
      PanelLimitersCs.Visible:=false;
   end;
end;

// ���������� ��� ������������ ������������� ��������� �������������.
procedure TFormSmagorinsky.CBLimitersCsClick(Sender: TObject);
begin
   // ��������� ��� ���������� ������������ ��������� �������������.
   if (CBLimitersCS.Checked) then
   begin
      // ������� Cs ��������
      Panel_user_limiters.Visible:=true;
      //EminCs.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].minCs);
      //EmaxCs.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].maxCs);
       EminCs.Text:=FloatToStr(Laplas.egddata.myflmod[0].minCs);
       EmaxCs.Text:=FloatToStr(Laplas.egddata.myflmod[0].maxCs);
   end
    else
   begin
      // ������� Cs ���������
      Panel_user_limiters.Visible:=false;
      //EminCs.Text:='-1e20';
      //EmaxCs.Text:='1e23';
   end;
end;

end.
