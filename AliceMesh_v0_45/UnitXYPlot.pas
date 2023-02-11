unit UnitXYPlot;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormXYPlot = class(TForm)
    GroupBoxXYPlot: TGroupBox;
    LabelXo: TLabel;
    LabelYo: TLabel;
    LabelZo: TLabel;
    Labeldirectional: TLabel;
    ButtonApply: TButton;
    EditXo: TEdit;
    EditYo: TEdit;
    EditZo: TEdit;
    ComboBoxdirectional: TComboBox;
    LabelXdim: TLabel;
    LabelYdim: TLabel;
    LabelZdim: TLabel;
    procedure ButtonApplyClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    bfirst_zapusk_XYPlot : Boolean;
  end;

var
  FormXYPlot: TFormXYPlot;

implementation

{$R *.dfm}

procedure TFormXYPlot.ButtonApplyClick(Sender: TObject);
begin
   Close;
end;

procedure TFormXYPlot.FormCreate(Sender: TObject);
begin
    bfirst_zapusk_XYPlot:=true; // Первый запуск.
end;

end.
