unit MeshUnit;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, Vcl.AppEvnts, Vcl.Buttons;

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
    BitBtnRegion: TBitBtn;
    procedure BApplyClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure BitBtnRegionClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  MeshForm: TMeshForm;

implementation
  uses
     VisualUnit, UnitAdaptRegionALICEMesh;
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

// Adaptive Region
procedure TMeshForm.BitBtnRegionClick(Sender: TObject);
begin
   case Laplas.ComboBoxlength.ItemIndex of
     0 : // m
     begin
        FormAdaptRedionALICEMesh.LabelUnitxS1.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitxE1.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnityS1.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnityE1.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitzS1.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitzE1.Caption:='m';

        FormAdaptRedionALICEMesh.LabelUnitxS2.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitxE2.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnityS2.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnityE2.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitzS2.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitzE2.Caption:='m';

        FormAdaptRedionALICEMesh.LabelUnitxS3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitxE3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnityS3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnityE3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitzS3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitzE3.Caption:='m';
     end;
     1 : // mm
     begin
        FormAdaptRedionALICEMesh.LabelUnitxS1.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitxE1.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnityS1.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnityE1.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitzS1.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitzE1.Caption:='mm';

        FormAdaptRedionALICEMesh.LabelUnitxS2.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitxE2.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnityS2.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnityE2.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitzS2.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitzE2.Caption:='mm';

        FormAdaptRedionALICEMesh.LabelUnitxS3.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitxE3.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnityS3.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnityE3.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitzS3.Caption:='mm';
        FormAdaptRedionALICEMesh.LabelUnitzE3.Caption:='mm';
     end;
     2 : // micron um
     begin
        FormAdaptRedionALICEMesh.LabelUnitxS1.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnitxE1.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnityS1.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnityE1.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnitzS1.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnitzE1.Caption:='um';

        FormAdaptRedionALICEMesh.LabelUnitxS2.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnitxE2.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnityS2.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnityE2.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnitzS2.Caption:='um';
        FormAdaptRedionALICEMesh.LabelUnitzE2.Caption:='um';

        FormAdaptRedionALICEMesh.LabelUnitxS3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitxE3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnityS3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnityE3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitzS3.Caption:='m';
        FormAdaptRedionALICEMesh.LabelUnitzE3.Caption:='m';
     end;
   end;
   FormAdaptRedionALICEMesh.ShowModal;
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
