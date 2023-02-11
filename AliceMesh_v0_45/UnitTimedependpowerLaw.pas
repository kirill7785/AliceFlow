unit UnitTimedependpowerLaw;
// ����������� �������� �������������� �� �������.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormTransientPowerSetting = class(TForm)
    Panel1: TPanel;
    RadioGroupTimeDependPowerLow: TRadioGroup;
    Button1: TButton;
    PanelTemperaturedependpower: TPanel;
    ComboBoxTemperaturedependpower: TComboBox;
    EditPower: TEdit;
    Label1: TLabel;
    Buttonpiecewisepower: TButton;
    procedure Button1Click(Sender: TObject);
    procedure FormClose(Sender: TObject; var Action: TCloseAction);
    procedure ButtonpiecewisepowerClick(Sender: TObject);
    procedure ComboBoxTemperaturedependpowerChange(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    bBlock_piecewice_power : Boolean;
  end;

var
  FormTransientPowerSetting: TFormTransientPowerSetting;

implementation

{$R *.dfm}
uses
     VisualUnit, addBlockUnit, UnitVariables, Unitusertempdepend;

// ���� ������
procedure TFormTransientPowerSetting.Button1Click(Sender: TObject);
var
  itype : Integer;
  // ��������������� ���������� ��� ��������� �������������� ��������
   bOk : Boolean;
   spow : String;
   rpow : Real;
begin
   // ����������� �������� �������������� �� �������.
   // 0 - �� ������� �� ������� � ���������� ���������,
   // 1 - squre wave ����������� �� �������,
   // 2 - square wave 2 ����������� �� �������,
   // 3 - hot cold ����� ��� ����������� �.�.
   // 4 - piecewise const.
   Laplas.body[Laplas.itek].ipower_time_depend:=RadioGroupTimeDependPowerLow.ItemIndex;
   bOk:=true; // ������� ������������ �����
   if (ComboBoxTemperatureDependPower.ItemIndex=0) then
   begin
      itype:=AddBlockForm.RadioGroupType.ItemIndex+1; // ��� �����
      case itype of
         1 : begin
                // SOLID
                spow:=EditPower.Text; // �������� ��������������
                if bOk then rpow:=FormVariables.my_real_convert(spow,bOk);
             end;
         2 : begin
                // HOLLOW
                spow:='0.0'; // �������� �� ����������
                rpow:=0.0;
             end;
         3 : begin
                // FLUID
                spow:=EditPower.Text; // �������� ��������������
                if bOk then rpow:=FormVariables.my_real_convert(spow,bOk);
             end;
       end; // case

       if (bOk) then
       begin
          Laplas.body[Laplas.itek].n_power:=1;
          SetLength(Laplas.body[Laplas.itek].temp_power, Laplas.body[Laplas.itek].n_power);
          SetLength(Laplas.body[Laplas.itek].arr_power, Laplas.body[Laplas.itek].n_power);
          SetLength(Laplas.body[Laplas.itek].arr_s_power, Laplas.body[Laplas.itek].n_power);
          Laplas.body[Laplas.itek].arr_s_power[0]:=spow;
          Laplas.body[Laplas.itek].arr_power[0]:=rpow;
          //AddBlockForm.labelpowerinfo.Caption:=FloatToStr(Laplas.body[Laplas.itek].arr_power[0])+' W';
          AddBlockForm.labelpowerinfo.Caption:=Trim((Laplas.body[Laplas.itek].arr_s_power[0])+' W');
          Close;
       end;
   end;
end;

// �������� �������������� ��������� �� �����������.
procedure TFormTransientPowerSetting.ButtonpiecewisepowerClick(Sender: TObject);
var
   i_4 : Integer;
begin
   if (not(bBlock_piecewice_power)) then
   begin
      Formusertempdepend.Caption:='Temperature depend power';
      Formusertempdepend.Label1.Caption:='The following curve specification consists of';
      Formusertempdepend.Label2.Caption:='a list of temperature/power pairs, which define a';
      Formusertempdepend.Label3.Caption:='piecewise-linear curve. Spacing is not significant';
      Formusertempdepend.Label4.Caption:='as long as the numbers are given in pairs.';
      Formusertempdepend.Label5.Caption:='Power units W.';
      Formusertempdepend.ComboBoxtemperatureUnit.ItemIndex:=0;  // ������� �������
      Formusertempdepend.Memopiecewiseproperties.Clear;
      for i_4 := 0 to Laplas.body[Laplas.itek].n_power-1 do
      begin
         Formusertempdepend.Memopiecewiseproperties.Lines.Add(FloatToStr(Laplas.body[Laplas.itek].temp_power[i_4])+' '+Laplas.body[Laplas.itek].arr_s_power[i_4]);
      end;
     Formusertempdepend.identifier:=2;   // power
     Formusertempdepend.ShowModal;
   end
   else
   begin
     // ���� �� ���� ��� ������ �� �������� �������� ���������
     // �� �����������.
      ComboBoxTemperaturedependpower.ItemIndex:=0;
      ComboBoxTemperaturedependpowerChange(Sender);
   end;
end;

procedure TFormTransientPowerSetting.ComboBoxTemperaturedependpowerChange(
  Sender: TObject);
begin
    if (ComboBoxTemperaturedependpower.ItemIndex=0) then
    begin
       // Constant
       FormTransientPowerSetting.EditPower.Visible:=true;
       FormTransientPowerSetting.Label1.Visible:=true;
       FormTransientPowerSetting.Buttonpiecewisepower.Visible:=false;
    end
     else
    begin
       // Piesewise Linear
       FormTransientPowerSetting.EditPower.Visible:=false;
       FormTransientPowerSetting.Label1.Visible:=false;
       FormTransientPowerSetting.Buttonpiecewisepower.Visible:=true;
    end;
end;



procedure TFormTransientPowerSetting.FormClose(Sender: TObject;
  var Action: TCloseAction);
begin
   //Laplas.body[Laplas.itek].ipower_time_depend:=RadioGroupTimeDependPowerLow.ItemIndex;
   Button1Click(Sender);
end;

procedure TFormTransientPowerSetting.FormCreate(Sender: TObject);
begin
    bBlock_piecewice_power:=false;
end;

end.
