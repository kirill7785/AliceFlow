unit Unitrectangularplot;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, TeEngine, Series, StdCtrls, ExtCtrls, TeeProcs, Chart,
  VclTee.TeeGDIPlus;

type
  TfrmRectangularPlot = class(TForm)
    cht1: TChart;
    lnsrsSeries1: TLineSeries;
    procedure FormResize(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  frmRectangularPlot: TfrmRectangularPlot;

implementation

{$R *.dfm}

// Изменение размеров экеранной формы.
procedure TfrmRectangularPlot.FormResize(Sender: TObject);
begin
   cht1.Height:=frmRectangularPlot.ClientHeight;
   cht1.Width:=frmRectangularPlot.ClientWidth;
end;

end.
