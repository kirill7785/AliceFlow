unit UnitSelProjMat;
// Здесь можно осуществить выбор уже существующего в проекте материала.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls;

type
  TFormSelProjMat = class(TForm)
    GroupBoxMain: TGroupBox;
    CBlistprojmat: TComboBox;
    BApply: TButton;
    procedure BApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormSelProjMat: TFormSelProjMat;

implementation
uses
     addBlockUnit, VisualUnit;
{$R *.dfm}

// по нажатию на кнопку осуществляется считывание выбранного материала
procedure TFormSelProjMat.BApplyClick(Sender: TObject);
begin
   Laplas.body[laplas.itek].imatid:=CBlistprojmat.ItemIndex;
   Close;
end;

end.
