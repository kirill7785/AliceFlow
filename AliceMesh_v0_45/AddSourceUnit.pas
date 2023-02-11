unit AddSourceUnit;
// ��������� ������� ��������� �����

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, ExtCtrls, StdCtrls;

type
  TAddSourceForm = class(TForm)
    Panelglobalcontainer: TPanel;
    RadioGroup1: TRadioGroup;
    Bapply: TButton;
    Panelinfo: TPanel;
    Lname: TLabel;
    Ename: TEdit;
    PanelGeometry: TPanel;
    GroupBox1: TGroupBox;
    LxS: TLabel;
    LyS: TLabel;
    LzS: TLabel;
    LxE: TLabel;
    LyE: TLabel;
    LzE: TLabel;
    ExS: TEdit;
    EyS: TEdit;
    EzS: TEdit;
    ExE: TEdit;
    EyE: TEdit;
    EzE: TEdit;
    RadioGroupPlane: TRadioGroup;
    PanelProperties: TPanel;
    GBpowerdef: TGroupBox;
    Label1: TLabel;
    LW: TLabel;
    RGpowertype: TRadioGroup;
    Ptempdefloc: TPanel;
    Label2: TLabel;
    Label3: TLabel;
    CBtableid: TComboBox;
    EOperoffsetdrain: TEdit;
    Epower: TEdit;
    PanelNetwork: TPanel;
    GroupBoxinx: TGroupBox;
    GroupBoxiny: TGroupBox;
    GroupBoxinz: TGroupBox;
    ComboBoxinx: TComboBox;
    ComboBoxiny: TComboBox;
    ComboBoxinz: TComboBox;
    procedure BapplyClick(Sender: TObject);
    procedure RadioGroupPlaneClick(Sender: TObject);
    procedure RGpowertypeClick(Sender: TObject);
    procedure RadioGroup1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  AddSourceForm: TAddSourceForm;

implementation
uses
     VisualUnit, UnitVariables;
{$R *.dfm}

procedure TAddSourceForm.BapplyClick(Sender: TObject);
var
   k : Integer;
   // ��������������� ���������� ��� ��������� �������������� ��������
   bOk : Boolean;
   s1, s2, s3, s4, s5, s6, spow : String;
   r1, r2, r3, r4, r5, r6, rpow : Real;
   buf : TPlane;

begin
   // ������������� :
   r1:=0.0;
   r2:=0.0;
   r3:=0.0;
   r4:=0.0;
   r5:=0.0;
   r6:=0.0;

   s1:=Trim(Epower.Text);
   (*
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
   *)
   if (FormatSettings.DecimalSeparator=',') then
      begin
         s1:=StringReplace(s1,'.',',',[rfReplaceAll]);
      end;
      if (FormatSettings.DecimalSeparator='.') then
      begin
         s1:=StringReplace(s1,',','.',[rfReplaceAll]);
      end;
   Epower.Text:=s1;


   s1:=Trim(ExS.Text);
   (*
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
   *)
    if (FormatSettings.DecimalSeparator='.') then
    begin
       s1:=StringReplace(s1,',','.',[rfReplaceAll]);
    end;

    if (FormatSettings.DecimalSeparator=',') then
    begin
       s1:=StringReplace(s1,'.',',',[rfReplaceAll]);
    end;
   ExS.Text:=s1;

    s1:=Trim(EyS.Text);
    (*
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
   *)
    if (FormatSettings.DecimalSeparator='.') then
    begin
      s1:=StringReplace(s1,',','.',[rfReplaceAll]);
    end;

    if (FormatSettings.DecimalSeparator=',') then
    begin
       s1:=StringReplace(s1,'.',',',[rfReplaceAll]);
    end;
   EyS.Text:=s1;

    s1:=Trim(EzS.Text);
    (*
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
   *)
    if (FormatSettings.DecimalSeparator='.') then
    begin
       s1:=StringReplace(s1,',','.',[rfReplaceAll]);
    end;

    if (FormatSettings.DecimalSeparator=',') then
    begin
       s1:=StringReplace(s1,'.',',',[rfReplaceAll]);
    end;
   EzS.Text:=s1;

   s1:=Trim(ExE.Text);
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
   ExE.Text:=s1;

    s1:=Trim(EyE.Text);
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
   EyE.Text:=s1;

    s1:=Trim(EzE.Text);
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
   EzE.Text:=s1;

    s1:=Trim(EOperoffsetdrain.Text);
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
   EOperoffsetdrain.Text:=s1;


   // ���� ������ �� ���������
   k:=Laplas.itek;
   //with (Laplas.source[k]) do
   buf:=Laplas.sourcepublic[k];
  //with (buf) do
   //begin
      bOk:=true; // ������� ������������ �����
      buf.iPlane:=RadioGroupPlane.ItemIndex+1;
      spow:=Epower.Text;
      if bOk then rpow:=FormVariables.my_real_convert(spow,bOk);
      // ��� ������� �������� : 0 - ���������, 1 - ��������.
      buf.itempdep:=RGpowertype.ItemIndex;
      if (buf.itempdep=1) then
      begin
         buf.soperatingoffsetdrain:=EOperoffsetdrain.Text;
         if bOk then buf.operatingoffsetdrain:=FormVariables.my_real_convert(buf.soperatingoffsetdrain,bOk);
         buf.id_table:=CBtableid.ItemIndex; // ���������� ����� �������
      end;

      buf.name:=Ename.Text; // ��� ��������
      // ������������� ����� ������� ����� �������� ����������� ���.
      Laplas.correctobjname('s',buf.name,k);
      Ename.Text:=buf.name;
      AddSourceForm.Caption:='Sources [ '+Trim(buf.name)+' ]';


      case buf.iPlane of
        1 : // XY
            begin
               s1:=ExS.Text;  // �����������������
               s2:=EyS.Text;  // ��������������
               s3:=ExE.Text;  // �������
               s4:=EyE.Text;  // ��������
               s5:=EzS.Text;  // �������������
               s6:=s5;

               if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // �������� �������
               if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // �������� �������������
               if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // � ������ �����������
               if bOk then r4:=FormVariables.my_real_convert(s4,bOk);  // �������� ����������.
               if bOk then r5:=FormVariables.my_real_convert(s5,bOk);
               if bOk then r6:=r5;
            end;
        2 : // XZ
            begin
               s1:=ExS.Text; // �����������������
               s2:=EyS.Text; // ��������������
               s3:=ExE.Text; // �������
               s4:=s2;      // ��������
               s5:=EzS.Text; // �������������
               s6:=EzE.Text;

               if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // �������� �������
               if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // �������� �������������
               if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // � ������ �����������
               if bOk then r4:=r2;  // �������� ����������.
               if bOk then r5:=FormVariables.my_real_convert(s5,bOk);
               if bOk then r6:=FormVariables.my_real_convert(s6,bOk);
            end;
        3 : // YZ
            begin
               s1:=ExS.Text; // �����������������
               s2:=EyS.Text; // ��������������
               s3:=s1;  // �������
               s4:=EyE.Text; // ��������
               s5:=EzS.Text; // �������������
               s6:=EzE.Text;

               if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // �������� �������
               if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // �������� �������������
               if bOk then r3:=r1; // � ������ �����������
               if bOk then r4:=FormVariables.my_real_convert(s4,bOk);  // �������� ����������.
               if bOk then r5:=FormVariables.my_real_convert(s5,bOk);
               if bOk then r6:=FormVariables.my_real_convert(s6,bOk);
            end;
      end; // case

      if (bOk) then
      begin
         buf.sxS:=s1;  // �����������������
         buf.syS:=s2;  // ��������������
         buf.sxE:=s3;  // �������
         buf.syE:=s4;  // ��������
         buf.szS:=s5;  // �������������
         buf.szE:=s6;
         buf.spower:=spow;

         buf.xS:=r1;  // �������� �������
         buf.yS:=r2;  // �������� �������������
         buf.xE:=r3;  // � ������ �����������
         buf.yE:=r4;  // �������� ����������.
         buf.zS:=r5;
         buf.zE:=r6;
         buf.Power:=rpow;

         buf.inx_network_conteiner:=ComboBoxinx.ItemIndex+1;
         buf.iny_network_conteiner:=ComboBoxiny.ItemIndex+1;
         buf.inz_network_conteiner:=ComboBoxinz.ItemIndex+1;
      end;

   //end;

   Laplas.sourcepublic[k]:=buf;
   with Laplas do
   begin
      ReadyPaint;
   end;

end;

procedure TAddSourceForm.RadioGroup1Click(Sender: TObject);
begin
   case RadioGroup1.ItemIndex of
     0 : begin
            // Info
            PanelInfo.Visible:=true;
            PanelGeometry.Visible:=false;
            PanelProperties.Visible:=false;
            PanelNetwork.Visible:=false;
         end;
     1 : begin
            // Geometry
            PanelInfo.Visible:=false;
            PanelGeometry.Visible:=true;
            PanelProperties.Visible:=false;
            PanelNetwork.Visible:=false;
         end;
     2 : begin
            // Properties
            PanelInfo.Visible:=false;
            PanelGeometry.Visible:=false;
            PanelProperties.Visible:=true;
            PanelNetwork.Visible:=false;
         end;
     3 : begin
            // Network
            PanelGeometry.Visible:=false;
            PanelInfo.Visible:=false;
            PanelProperties.Visible:=false;
            GroupBoxinx.Visible:=true;
            GroupBoxiny.Visible:=true;
            GroupBoxinz.Visible:=true;
            case RadioGroupPlane.ItemIndex of
              0 : begin
                     // XY
                    GroupBoxinz.Visible:=false;
                    ComboBoxinz.ItemIndex:=0;
              end;
              1 : begin
                     // XZ
                    GroupBoxiny.Visible:=false;
                    ComboBoxiny.ItemIndex:=0;
              end;
              2 : begin
                     // YZ
                    GroupBoxinx.Visible:=false;
                    ComboBoxinx.ItemIndex:=0;
              end;
            end;
            PanelNetwork.Visible:=true;
         end;
   end;
end;

procedure TAddSourceForm.RadioGroupPlaneClick(Sender: TObject);
begin
   // ����� ��������� � ������� ����� �������� �����
   case (RadioGroupPlane.ItemIndex+1) of
     1 : // XY
         begin
            ExE.Visible:=true;
            LxE.Visible:=true;
            EyE.Visible:=true;
            LyE.Visible:=true;
            EzE.Visible:=false;
            LzE.Visible:=false;
         end;
     2 : // XZ
         begin
            ExE.Visible:=true;
            LxE.Visible:=true;
            EyE.Visible:=false;
            LyE.Visible:=false;
            EzE.Visible:=true;
            LzE.Visible:=true;
         end;
     3 : // YZ
         begin
            ExE.Visible:=false;
            LxE.Visible:=false;
            EyE.Visible:=true;
            LyE.Visible:=true;
            EzE.Visible:=true;
            LzE.Visible:=true;
         end;
   end;
end;

procedure TAddSourceForm.RGpowertypeClick(Sender: TObject);
var
    i : Integer;
begin
    // ����� ������� ������� �������� :
    // �������� ����� ���������� ���� ���������� (���������)
    // ���� � ���� ������� ��� ����������� �� ����������� � �������� �����.
    case RGpowertype.ItemIndex of
       0 : // const
           begin
              Ptempdefloc.Visible:=false;
              Label1.Caption:='power';
              LW.Caption:='W';
           end;
       1 : // power define
           begin
              if (Laplas.iltdp=0) then
              begin
                  Application.MessageBox('Plese Define -> Power Table create','Please define',MB_OK);
                  Laplas.MainMemo.Lines.Add('Plese Define -> Power Table create');
                  RGpowertype.ItemIndex:=0; // ���������� ��������.
                  Ptempdefloc.Visible:=false;
                  Label1.Caption:='power';
                  LW.Caption:='W';
              end
               else
              begin
                 // ���������� ������ ��������� ������.
                 CBtableid.Clear;
                 for i:=0 to Laplas.iltdp-1 do
                 begin
                    CBtableid.AddItem(IntToStr(i),Sender);
                 end;
                 CBtableid.ItemIndex:=Laplas.source[Laplas.itek].id_table;  // ���������� ����� �������� �������� ��������.
                 EOperoffsetdrain.Text:=Laplas.source[Laplas.itek].soperatingoffsetdrain; // ���������� �� �����
                 Ptempdefloc.Visible:=true;
                 // mult power - ��� ���������
                 // �� ������� ����������� �������� ��������.
                 // ��� ������ ���� ��������� ������ �������� ��������� ���������
                 // � ������������ ����������, ����� �������� ���� ��������� �������
                 // ��� ������������� mult power ������ 0.5.
                 Label1.Caption:='mult power';
                 LW.Caption:='';
              end;
           end;
    end;
end;

end.

