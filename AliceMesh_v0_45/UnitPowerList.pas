unit UnitPowerList;
// ���� ������ ������������ ��� ������������
// ������ ��� ������. � ������ ���������� �����
// � ���������� ������ ���������� ������� �����
// � ������� �������� ������ ����������� ������������
// �������� �� ����������� � �������� �����.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TFormPowerList = class(TForm)
    Panelmain: TPanel;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Eltdp: TEdit;
    BApplyglobal: TButton;
    GBtableparam: TGroupBox;
    Label5: TLabel;
    CBidtable: TComboBox;
    Label6: TLabel;
    Efilename: TEdit;
    Bonetable: TButton;
    procedure BApplyglobalClick(Sender: TObject);
    procedure BonetableClick(Sender: TObject);
    procedure CBidtableChange(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormPowerList: TFormPowerList;

implementation
 uses
     VisualUnit;

{$R *.dfm}

procedure TFormPowerList.BApplyglobalClick(Sender: TObject);
var
    i, innew : Integer;
begin
   // ������� ���������� ������

   innew:=StrToInt(Eltdp.Text);
   if (innew <> Laplas.iltdp) then
   begin
       // ����� ������ � ������ ��������� �������.
       Laplas.iltdp:=innew;
       SetLength(Laplas.listtablename, Laplas.iltdp);
       CBidtable.Clear;
       for i:=0 to Laplas.iltdp-1 do
       begin
          CBidtable.AddItem(IntToStr(i),Sender);
          Laplas.listtablename[i]:='';
       end;
       if (Laplas.iltdp>0) then
       begin
          CBidtable.ItemIndex:=0;
          GBtableparam.Visible:=true;
       end
       else GBtableparam.Visible:=false;
       Efilename.Text:=''; // ������ ������.
   end;
end;

procedure TFormPowerList.BonetableClick(Sender: TObject);
begin
   // ������� ����� ���������� �������:
   Laplas.listtablename[CBidtable.ItemIndex]:=Efilename.Text;
end;

procedure TFormPowerList.CBidtableChange(Sender: TObject);
begin
   // ����� ������ �������:
   Efilename.Text:=Laplas.listtablename[CBidtable.ItemIndex];
end;

end.
