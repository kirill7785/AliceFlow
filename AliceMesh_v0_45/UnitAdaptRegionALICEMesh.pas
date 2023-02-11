unit UnitAdaptRegionALICEMesh;
// Адаптация региона  в форме прямоугольного параллепипеда для АЛИС сетки.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormAdaptRedionALICEMesh = class(TForm)
    GroupBoxCabinet: TGroupBox;
    ComboBoxCabinetMinLevel: TComboBox;
    Label1: TLabel;
    GroupBoxAdaptiveRegion1: TGroupBox;
    LabelNameminLevel1: TLabel;
    ComboBoxAdaptiveRegion1: TComboBox;
    LabelNamexS1: TLabel;
    EditxS1: TEdit;
    LabelUnitxS1: TLabel;
    LabelNamexE1: TLabel;
    EditxE1: TEdit;
    LabelUnitxE1: TLabel;
    LabelNameyS1: TLabel;
    EdityS1: TEdit;
    LabelUnityS1: TLabel;
    LabelNameyE1: TLabel;
    EdityE1: TEdit;
    LabelUnityE1: TLabel;
    LabelNamezS1: TLabel;
    EditzS1: TEdit;
    LabelUnitzS1: TLabel;
    LabelNamezE1: TLabel;
    EditzE1: TEdit;
    LabelUnitzE1: TLabel;
    GroupBoxAdaptiveRegion2: TGroupBox;
    LabelNameminLevel2: TLabel;
    LabelNamexS2: TLabel;
    LabelUnitxS2: TLabel;
    LabelNamexE2: TLabel;
    LabelUnitxE2: TLabel;
    LabelNameyS2: TLabel;
    LabelUnityS2: TLabel;
    LabelNameyE2: TLabel;
    LabelUnityE2: TLabel;
    LabelNamezS2: TLabel;
    LabelUnitzS2: TLabel;
    LabelNamezE2: TLabel;
    LabelUnitzE2: TLabel;
    ComboBoxAdaptiveRegion2: TComboBox;
    EditxS2: TEdit;
    EditxE2: TEdit;
    EdityS2: TEdit;
    EdityE2: TEdit;
    EditzS2: TEdit;
    EditzE2: TEdit;
    GroupBoxAdaptiveRegion3: TGroupBox;
    LabelNameminLevel3: TLabel;
    LabelNamexS3: TLabel;
    LabelUnitxS3: TLabel;
    LabelNamexE3: TLabel;
    LabelUnitxE3: TLabel;
    LabelNameyS3: TLabel;
    LabelUnityS3: TLabel;
    LabelNameyE3: TLabel;
    LabelUnityE3: TLabel;
    LabelNamezS3: TLabel;
    LabelUnitzS3: TLabel;
    LabelNamezE3: TLabel;
    LabelUnitzE3: TLabel;
    ComboBoxAdaptiveRegion3: TComboBox;
    EditxS3: TEdit;
    EditxE3: TEdit;
    EdityS3: TEdit;
    EdityE3: TEdit;
    EditzS3: TEdit;
    EditzE3: TEdit;
    ButtonApply: TButton;
    procedure ButtonApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormAdaptRedionALICEMesh: TFormAdaptRedionALICEMesh;

implementation

{$R *.dfm}

uses VisualUnit;

// Адаптация региона  в форме прямоугольного параллепипеда для АЛИС сетки.
procedure TFormAdaptRedionALICEMesh.ButtonApplyClick(Sender: TObject);
begin

   if (ComboBoxAdaptiveRegion3.ItemIndex>0) then
   begin
       Laplas.isizearrAdaptRegion:=3;
   end
   else if (ComboBoxAdaptiveRegion2.ItemIndex>0) then
   begin
      Laplas.isizearrAdaptRegion:=2;
   end
   else if (ComboBoxAdaptiveRegion1.ItemIndex>0) then
   begin
      Laplas.isizearrAdaptRegion:=1;
   end
   else
   begin
      Laplas.isizearrAdaptRegion:=0;
   end;

   Laplas.ilevelMeshCabinet:=ComboBoxCabinetMinLevel.ItemIndex-1;

   Close;
end;

end.
