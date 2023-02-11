object FormRadiation: TFormRadiation
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Radiation inside the block'
  ClientHeight = 122
  ClientWidth = 467
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBoxemissivity: TGroupBox
    Left = 0
    Top = 0
    Width = 467
    Height = 122
    Caption = 'emissivity individual side of block'
    TabOrder = 0
    object Label6: TLabel
      Left = 16
      Top = 24
      Width = 82
      Height = 13
      Caption = 'epsilon min X (W)'
    end
    object Label7: TLabel
      Left = 16
      Top = 56
      Width = 82
      Height = 13
      Caption = 'epsilon max X (E)'
    end
    object Label8: TLabel
      Left = 151
      Top = 24
      Width = 78
      Height = 13
      Caption = 'epsilon min Y (S)'
    end
    object Label9: TLabel
      Left = 151
      Top = 56
      Width = 83
      Height = 13
      Caption = 'epsilon max Y (N)'
    end
    object Label10: TLabel
      Left = 297
      Top = 24
      Width = 78
      Height = 13
      Caption = 'epsilon min Z (B)'
    end
    object Label11: TLabel
      Left = 293
      Top = 64
      Width = 82
      Height = 13
      Caption = 'epsilon max Z (T)'
    end
    object Edit1: TEdit
      Left = 104
      Top = 16
      Width = 41
      Height = 21
      TabOrder = 0
      Text = '0.8'
    end
    object Edit2: TEdit
      Left = 104
      Top = 58
      Width = 41
      Height = 21
      TabOrder = 1
      Text = '0.8'
    end
    object Edit3: TEdit
      Left = 235
      Top = 16
      Width = 46
      Height = 21
      TabOrder = 2
      Text = '0.8'
    end
    object Edit4: TEdit
      Left = 235
      Top = 56
      Width = 52
      Height = 21
      TabOrder = 3
      Text = '0.8'
    end
    object Edit5: TEdit
      Left = 381
      Top = 16
      Width = 67
      Height = 21
      TabOrder = 4
      Text = '0.8'
    end
    object Edit6: TEdit
      Left = 381
      Top = 56
      Width = 67
      Height = 21
      TabOrder = 5
      Text = '0.8'
    end
    object ButtonApply: TButton
      Left = 373
      Top = 83
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 6
      OnClick = ButtonApplyClick
    end
    object CheckBoxinternalRadiation: TCheckBox
      Left = 24
      Top = 96
      Width = 105
      Height = 17
      Caption = 'Internal Radiation'
      TabOrder = 7
    end
  end
end
