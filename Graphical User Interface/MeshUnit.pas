unit MeshUnit;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, Vcl.AppEvnts;

type
  TMeshForm = class(TForm)
    GBnumnodes: TGroupBox;
    Linx: TLabel;
    Liny: TLabel;
    Linz: TLabel;
    BApply: TButton;
    CBinx: TComboBox;
    CBiny: TComboBox;
    CBinz: TComboBox;
    lbl1: TLabel;
    edtmaxsizeratio: TEdit;
    ComboBoxmeshgen: TComboBox;
    CheckBoxALICE: TCheckBox;
    Label1: TLabel;
    ApplicationEvents1: TApplicationEvents;
    Editmaxsizeratio2: TEdit;
    ComboBoxALICEType: TComboBox;
    ComboBoxSnapTo: TComboBox;
    procedure BApplyClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  MeshForm: TMeshForm;

implementation
  uses
     VisualUnit;
{$R *.dfm}

// Запрет форме сворачиваться.
procedure TMeshForm.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
     if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
   then
      msg.message:=0;
end;

// Ввод параметров при нажатии кнопки Apply-применить.
procedure TMeshForm.BApplyClick(Sender: TObject);
var
  s : string;
  i : Integer;
begin
   Laplas.inx:=CBinx.ItemIndex+4;
   Laplas.iny:=CBiny.ItemIndex+4;
   Laplas.inz:=CBinz.ItemIndex+4;
   s:=Trim(edtmaxsizeratio.Text);
   if (FormatSettings.DecimalSeparator='.') then
   begin
      for i:=1 to Length(s) do
      begin
         if (s[i]=',') then
         begin
            s[i]:='.';
         end;
      end;
   end;
   if (FormatSettings.DecimalSeparator=',') then
   begin
      for i:=1 to Length(s) do
      begin
         if (s[i]='.') then
         begin
            s[i]:=',';
         end;
      end;
   end;
   Laplas.etalon_max_size_ratio:=StrToFloat(s);
   if ( Laplas.etalon_max_size_ratio<=1.0) then
   begin
      Laplas.etalon_max_size_ratio:=2.0;
      edtmaxsizeratio.Text:='2.0';
   end;
   s:=Trim(Editmaxsizeratio2.Text);
   if (FormatSettings.DecimalSeparator='.') then
   begin
      for i:=1 to Length(s) do
      begin
         if (s[i]=',') then
         begin
            s[i]:='.';
         end;
      end;
   end;
   if (FormatSettings.DecimalSeparator=',') then
   begin
      for i:=1 to Length(s) do
      begin
         if (s[i]='.') then
         begin
            s[i]:=',';
         end;
      end;
   end;
   Laplas.etalon_max_size_ratio2:=StrToFloat(s);
   Close();
end;

// Инициализация при создании формы
procedure TMeshForm.FormCreate(Sender: TObject);
begin
   CBinx.ItemIndex:=Laplas.inx-4;
   CBiny.ItemIndex:=Laplas.iny-4;
   CBinz.ItemIndex:=Laplas.inz-4;
   edtmaxsizeratio.Text:=FloatToStr(Laplas.etalon_max_size_ratio);
end;

end.
