object FormVariables: TFormVariables
  Left = 81
  Top = 170
  AutoSize = True
  Caption = 'Add or redo variables'
  ClientHeight = 377
  ClientWidth = 433
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object GBVariables: TGroupBox
    Left = 0
    Top = 0
    Width = 433
    Height = 377
    Caption = 'Add variables to project'
    Color = clMoneyGreen
    ParentColor = False
    TabOrder = 0
    object StringGridVariables: TStringGrid
      Left = 8
      Top = 24
      Width = 393
      Height = 305
      ColCount = 3
      DefaultColWidth = 120
      RowCount = 12
      Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goDrawFocusSelected, goRowSizing, goColSizing, goEditing, goTabs]
      TabOrder = 0
    end
    object BApply: TButton
      Left = 304
      Top = 335
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = BApplyClick
    end
    object ButtonRename: TButton
      Left = 152
      Top = 335
      Width = 75
      Height = 25
      Caption = 'Rename'
      TabOrder = 2
      OnClick = ButtonRenameClick
    end
  end
end
